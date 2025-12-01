create table guiders (
	guiders_id       	 bigserial primary key,    
	dateobs          varchar ,
	acquired          boolean    ,
	guiding          boolean    ,
	guide_target_x		 float    ,
	guide_target_y		 float    ,
	exp_time	 	 float    ,
	exp_avg   	 	 float    ,
	ins_at                  timestamp  not null
);

