create table weather (
	weather_station_id    bigserial primary key,    
	dateobs                varchar, 
        temp                   real,
        hum                    real,
        pressure               real,
        dewpoint               real,
        wind_speed             real,
        wind_dir               real,
        encl35m                integer,
        encl25m                integer
);
