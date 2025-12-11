create table motors (
	motors_id       bigserial primary key,    
	iodine_pos           	 real,
	calib_mirror_pos	 real,
	spectrograph_foc	 int    ,
	iodine1_temp_set	 	 real    ,
	iodine1_temp_read	 real    ,
	iodine1_current 	 real    ,
	iodine1_voltage 	 real    ,
	iodine2_temp_set	 	 real    ,
	iodine2_temp_read	 real    ,
	iodine2_current 	 real    ,
	iodine2_voltage 	 real    ,
	lamp_halogen_on	 	 int    ,
	lamp_thar_on	 	 int    ,
	lamp_led_on	         int,
	ins_at                   timestamp  not null
);

