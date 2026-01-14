CREATE SEQUENCE guiders_guider_id_seq;

CREATE TABLE IF NOT EXISTS public.guiders
(
    guiders_id bigint NOT NULL DEFAULT nextval('guiders_guiders_id_seq'::regclass),
    dateobs character varying COLLATE pg_catalog."default",
    acquired boolean,
    guiding boolean,
    guide_target_x double precision,
    guide_target_y double precision,
    exp_time double precision,
    exp_avg double precision,
    ins_at timestamp without time zone NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT guiders_pkey PRIMARY KEY (guiders_id)
)

TABLESPACE pg_default;

ALTER TABLE IF EXISTS public.guiders
    OWNER to song;


-- Table: public.motors

-- DROP TABLE IF EXISTS public.motors;

CREATE SEQUENCE motors_motors_id_seq;
CREATE TABLE IF NOT EXISTS public.motors
(
    motors_id bigint NOT NULL DEFAULT nextval('motors_motors_id_seq'::regclass),
    iodine_pos real,
    calib_mirror_pos real,
    spectrograph_foc integer,
    iodine1_temp_set real,
    iodine1_temp_read real,
    iodine1_current real,
    iodine1_voltage real,
    iodine2_temp_set real,
    iodine2_temp_read real,
    iodine2_current real,
    iodine2_voltage real,
    lamp_halogen_on integer,
    lamp_thar_on integer,
    lamp_led_on integer,
    ins_at timestamp without time zone NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT motors_pkey PRIMARY KEY (motors_id)
)

TABLESPACE pg_default;

ALTER TABLE IF EXISTS public.motors
    OWNER to song;

-- Table: public.tel_dome

-- DROP TABLE IF EXISTS public.tel_dome;
CREATE SEQUENCE tel_dome_tel_dome_id_seq;
CREATE TABLE IF NOT EXISTS public.tel_dome
(
    tel_dome_id bigint NOT NULL DEFAULT nextval('tel_dome_tel_dome_id_seq'::regclass),
    tel_ready_state integer,
    tel_con_state boolean,
    tel_tracking boolean,
    tel_ra_j2000 real,
    tel_dec_j2000 real,
    tel_ra real,
    tel_dec real,
    tel_alt real,
    tel_azm real,
    tel_alt_rms_error real,
    tel_azm_rms_error real,
    m3_pos integer,
    dome_shutterstate integer,
    dome_az real,
    dome_slewing boolean,
    dome_light_state boolean,
    temp_m1 real,
    temp_m2 real,
    temp_m3 real,
    temp_back real,
    temp_amb real,
    focuser_1_pos integer,
    focuser_1_moving boolean,
    focuser_2_pos integer,
    focuser_2_moving boolean,
    tel_lst character varying COLLATE pg_catalog."default",
    ins_at timestamp without time zone NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT tel_dome_pkey PRIMARY KEY (tel_dome_id)
)

TABLESPACE pg_default;

ALTER TABLE IF EXISTS public.tel_dome
    OWNER to song;

-- Table: public.weather

-- DROP TABLE IF EXISTS public.weather;

CREATE SEQUENCE weather_weather_station_id_seq;
CREATE TABLE IF NOT EXISTS public.weather
(
    weather_id bigint NOT NULL DEFAULT nextval('weather_weather_station_id_seq'::regclass),
    dateobs character varying COLLATE pg_catalog."default",
    temp real,
    hum real,
    pressure real,
    dewpoint real,
    wind_speed real,
    wind_dir real,
    encl35m integer,
    encl25m integer,
    ins_at timestamp with time zone NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT weather_pkey PRIMARY KEY (weather_id)
)

TABLESPACE pg_default;

ALTER TABLE IF EXISTS public.weather
    OWNER to song;
