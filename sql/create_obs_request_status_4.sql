create table obs_request_status4 (
        req_no          serial primary key,
        status		 varchar not null,
        no_exp		 int,
        comment	 varchar,
        ins_at          timestamp not null default current_timestamp
);
