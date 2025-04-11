import os
import influxdb_client
from influxdb_client.client.write_api import SYNCHRONOUS
import pdb

os.environ['INFLUX_TOKEN'] = 'v-RuHY6T1pyOIL1SU9lrWYKYEU_SDZ0VWkPHOIU9hMECF7axu2wiFzY1u8N7J6s9KCbOreQKI43mJUi9pj5BbA=='

os.environ['INFLUX_TOKEN'] = 'unP6W2cm9uLUWIi9jVj8tDjhYsN286L3NWcmuNFN8tvJOh9Pewgc9LQLPYx302XAebNhHIVnv0ubi365V0oUPA=='


# Store the URL of your InfluxDB instance
url="http://localhost:8086"
token = os.environ['INFLUX_TOKEN']
org='NMSU'
client = influxdb_client.InfluxDBClient(
            url=url,
            token=token,
            org=org
         )
write_api = client.write_api(write_options=SYNCHRONOUS)

def write(idict,bucket=None,measurement=None,org='NMSU',location='apo') :
    """ Add measurements to influx database to desired bucket,measurement
    """
    p = []
    for key in idict.keys() :
        p.append(influxdb_client.Point(measurement).tag("location",location).field(key,idict[key]))
    write_api.write(bucket=bucket, org=org, record=p)

