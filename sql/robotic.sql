/*GRANT CREATE ON DATABASE aposong To song ;*/
CREATE TABLE robotic.target 
(
    target_pk bigserial primary key, 
    targname text,
    ra text,
    dec text,
    epoch double precision
);

CREATE TABLE robotic.schedule (
   schedule_pk bigserial primary key,
   schedulename text,
   min_airmass double precision,
   max_airmass double precision,
   nvisits integer,
   dt_visit double precision,
   nsequence integer
) ;

CREATE TABLE robotic.sequence (
   sequence_pk bigserial primary key,
   sequencename text,
   n_exp integer[],
   t_exp double precision[],
   filter text[]
) ;


CREATE TABLE robotic.focus (
   focus_pk bigserial primary key,
   start_foc integer,
   end_foc integer,
   date integer,
   mjd double precision,
   azimuth double precision,
   altitude double precision,
   best double precision,
   best_fit double precision 
) ;
)
