CREATE TABLE tags (
id     serial PRIMARY KEY,
tag    text,
UNIQUE(tag)
);

CREATE TABLE annotations (
id     serial PRIMARY KEY,
tags   integer[]
);

CREATE TABLE kmers (
id              serial PRIMARY KEY,
hash            integer,
annotation_id   integer,
FOREIGN KEY (annotation_id) REFERENCES annotations (id)
);
