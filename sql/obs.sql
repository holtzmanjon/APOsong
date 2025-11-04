/*GRANT CREATE ON DATABASE aposong To song ;*/
CREATE TABLE obs.exposure
(
    exp_pk bigserial primary key, 
    dateobs text,
    mjd double precision,
    ra double precision,
    dec double precision,
    az double precision,
    alt double precision,
    rot double precision,
    focus double precision,
    exptime double precision,
    filter text,
    ccdtemp double precision,
    xbin integer,
    ybin integer,
    file text
);

CREATE TABLE obs.focus (
    focus_pk bigserial primary key,
    mjd double precision,
    exptime double precision,
    filter text,
    camera integer,
    bin integer,
    focvals integer[],
    bestfoc double precision,
    besthf double precision,
    bestfitfoc double precision,
    bestfithf double precision,
    files text[]
);

CREATE TABLE obs.reduced (
    reduced_pk bigserial primary key,
    exp_pk integer,
    throughput double precision,
    sn double precision
);
