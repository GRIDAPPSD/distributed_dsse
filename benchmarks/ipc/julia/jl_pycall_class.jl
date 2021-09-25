#!/usr/bin/env julia

using PyCall

py"""
from gridappsd import GridAPPSD

hello = 'Sai'

class AMQWrapper:
  def __init__(self):
    self.gapps = GridAPPSD()
    self.reply = []

  def bye(self):
    return hello

  def query(self, querystr):
    results = self.gapps.query_data(querystr)
    return results['data']['results']['bindings']


  def callback(self, header, message):
    #print('Python callback topic: ' + header['destination'])
    if 'message' in message:
      self.reply.append((header['destination'], message['message']))
    else:
      self.reply.append((header['destination'], message))


  def subscribe(self, topic):
    self.gapps.subscribe(topic, self.callback)


  def clearMessage(self):
    self.reply.clear()
    return self.reply


  def getMessage(self):
    #print('Python reply count: ' + str(len(self.reply)))
    return self.reply
"""

CallbackDict = Dict()

function logCallback(message)
  println(string("Julia log status: ", message["processStatus"]))
  completeFlag = message["processStatus"]=="COMPLETE"
  return !completeFlag
end


function measCallback(message)
  println(string("Julia measurement timestamp: ", message["timestamp"]))
  for (ment_mrid, ment_data) in message["measurements"]
    mag = "Unspecified"
    if haskey(ment_data, "magnitude")
      mag = ment_data["magnitude"]
    end
    ang = "Unspecified"
    if haskey(ment_data, "angle")
      ang = ment_data["angle"]
    end
    println(string("Julia measurement mrid: ", ment_mrid, ", magnitude: ", mag, ", angle: ", ang))
  end
  return true
end


function AMQSubscribe(pyamq, topic, callback)
  CallbackDict[topic] = callback
  pyamq.subscribe(topic)
end


function AMQMainLoop(pyamq)
  keepLoopingFlag = true

  while keepLoopingFlag
    reply = pyamq.clearMessage()

    while length(reply) == 0
      sleep(0.1)
      reply = pyamq.getMessage()
    end

    for (topic, message) in reply
      if !CallbackDict[topic](message)
        keepLoopingFlag = false
      end
    end
  end
end


# main

pyamq = py"AMQWrapper"()

feeder_mrid = "_5B816B93-7A5F-B64C-8460-47C17D6E4B0F"

SOURCE_QUERY =
  "PREFIX r:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n"*
  "PREFIX c:  <http://iec.ch/TC57/CIM100#>\n"*
  "SELECT DISTINCT ?busname ?nomv WHERE {\n"*
  "VALUES ?fdrid {\""*
  feeder_mrid*
  "\"}\n"*
  "?fdr c:IdentifiedObject.mRID ?fdrid.\n"*
  "?bus c:ConnectivityNode.ConnectivityNodeContainer ?fdr.\n"*
  "?bus r:type c:ConnectivityNode.\n"*
  "?bus c:IdentifiedObject.name ?busname.\n"*
  "?bus c:IdentifiedObject.mRID ?cnid.\n"*
  "?fdr c:IdentifiedObject.name ?feeder.\n"*
  "?trm c:Terminal.ConnectivityNode ?bus.\n"*
  "?trm c:Terminal.ConductingEquipment ?ce.\n"*
  "?ce c:ConductingEquipment.BaseVoltage ?bv.\n"*
  "?bv c:BaseVoltage.nominalVoltage ?nomv.\n"*
  "}\n"*
  "ORDER by ?feeder ?busname ?nomv"

results = pyamq.query(SOURCE_QUERY)
for obj in results
  println(string("Julia busname: ", obj["busname"]["value"], ", nomv: ", obj["nomv"]["value"]))
end

what = pyamq.bye()
println(string("What: ", what))


if length(ARGS) > 0
  sim_id = ARGS[1]

  AMQSubscribe(pyamq, "/topic/goss.gridappsd.simulation.log."*sim_id, logCallback)
  AMQSubscribe(pyamq, "/topic/goss.gridappsd.simulation.output."*sim_id, measCallback)

  AMQMainLoop(pyamq)
end

println("Julia all done!")

