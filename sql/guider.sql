create table guiders (
	guiders_id       	 bigserial primary key,    
	guiding_active          boolean    ,
	paused   	        boolean    ,
	pointing_enabled	 boolean ,
	guide_target_x		 float    ,
	guide_target_y		 float    ,
	exp_time	 	 float    ,
	slow_sig		 float,
	fast_sig		 float,
	found_star		 boolean    ,
	star_on_fibre		 boolean,
	seeing			 float,
	flux			 float,
	number_of_stars	 int,
	extra_val_1		 real,
	extra_val_2		 real,
	extra_val_3		 int,
	extra_val_4		 int,
	extra_param_1		 varchar,
	extra_param_2		 varchar,
	ins_at                  timestamp  not null
);

