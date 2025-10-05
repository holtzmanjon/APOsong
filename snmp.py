from pysnmp.carrier.asyncio.dispatch import AsyncioDispatcher
from pysnmp.carrier.asyncio.dgram import udp, udp6
from pyasn1.codec.ber import encoder, decoder
from pysnmp.proto import api
from pysnmp.proto.rfc1905 import ResponsePDU
import pdb
import time
import influx
import numpy as np

out = None

def snmpget(oid,host) :
    global out

    # Protocol version to use
    pMod = api.PROTOCOL_MODULES[api.SNMP_VERSION_1]
    # pMod = api.PROTOCOL_MODULES[api.SNMP_VERSION_2C]
    # Build PDU
    reqPDU = pMod.GetRequestPDU()
    pMod.apiPDU.set_defaults(reqPDU)
    pMod.apiPDU.set_varbinds(
            reqPDU, ((oid, pMod.Null("")),
                     (oid, pMod.Null("")))
            )

    # Build message
    reqMsg = pMod.Message()
    pMod.apiMessage.set_defaults(reqMsg)
    pMod.apiMessage.set_community(reqMsg, "public")
    pMod.apiMessage.set_pdu(reqMsg, reqPDU)


    # noinspection PyUnusedLocal,PyUnusedLocal
    def cbRecvFun(
        transportDispatcher, transportDomain, transportAddress, wholeMsg, reqPDU=reqPDU
    ):
        global out
        while wholeMsg:
            rspMsg, wholeMsg = decoder.decode(wholeMsg, asn1Spec=pMod.Message())
            rspPDU: ResponsePDU = pMod.apiMessage.get_pdu(rspMsg)

            # Match response to request
            if pMod.apiPDU.get_request_id(reqPDU) == pMod.apiPDU.get_request_id(rspPDU):
                # Check for SNMP errors reported
                errorStatus = pMod.apiPDU.get_error_status(rspPDU)
                if errorStatus:
                    print(errorStatus.prettyPrint())

                else:
                    for oid, val in pMod.apiPDU.get_varbinds(rspPDU):
                        print(f"{oid.prettyPrint()} = {val.prettyPrint()}")
                        out = val.prettyPrint()

                transportDispatcher.job_finished(1)

        return wholeMsg


    transportDispatcher = AsyncioDispatcher()

    transportDispatcher.register_recv_callback(cbRecvFun)

    # UDP/IPv4
    transportDispatcher.register_transport(
        udp.DOMAIN_NAME, udp.UdpAsyncioTransport().open_client_mode()
    )

    # Pass message to dispatcher
    transportDispatcher.send_message(
        encoder.encode(reqMsg), udp.DOMAIN_NAME, (host, 161)
    )
    transportDispatcher.job_started(1)

## UDP/IPv6 (second copy of the same PDU will be sent)
# transportDispatcher.register_transport(
#    udp6.domainName, udp6.Udp6AsyncioTransport().open_client_mode()
# )

# Pass message to dispatcher
# transportDispatcher.send_message(
#    encoder.encode(reqMsg), udp6.domainName, ('::1', 161)
# )
# transportDispatcher.job_started(1)

# Dispatcher will finish as job#1 counter reaches zero
    transportDispatcher.run_dispatcher(3)

    transportDispatcher.close_dispatcher()

def ftoc(x) :
    return (int(x)/100.-32.)*5/9

offsets=np.zeros([3,2])
offsets[0,0]= 0.972 - 0.06
offsets[1,0]= 0.972-1.367 - 0.06
offsets[2,0]= 0. - 0.06
offsets[0,1]= 0.967 + 0.06
offsets[1,1]= 0. + 0.06
offsets[2,1]= 0.967-1.505 + 0.06
while True :
    dict={}
    for i,host in enumerate(["10.75.0.18","10.75.0.25"]) :
        try :
            snmpget("1.3.6.1.4.1.20916.1.7.1.1.1.2.0",host)
            dict['t1']=ftoc(out)+offsets[0,i]
            snmpget("1.3.6.1.4.1.20916.1.7.1.2.1.2.0",host)
            dict['t2']=ftoc(out)+offsets[1,i]
            snmpget("1.3.6.1.4.1.20916.1.7.1.2.2.2.0",host)
            dict['t3']=ftoc(out)+offsets[2,i]

            influx.write(dict,bucket='spectemp',measurement=f'tempager_{i}')
        except : pass

    time.sleep(60)
