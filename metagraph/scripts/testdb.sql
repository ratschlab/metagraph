--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET search_path = public, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: annotations; Type: TABLE; Schema: public; Owner: gideon; Tablespace: 
--

CREATE TABLE annotations (
    id integer NOT NULL,
    tags integer[]
);


ALTER TABLE annotations OWNER TO gideon;

--
-- Name: annotations_id_seq; Type: SEQUENCE; Schema: public; Owner: gideon
--

CREATE SEQUENCE annotations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE annotations_id_seq OWNER TO gideon;

--
-- Name: annotations_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gideon
--

ALTER SEQUENCE annotations_id_seq OWNED BY annotations.id;


--
-- Name: kmers; Type: TABLE; Schema: public; Owner: gideon; Tablespace: 
--

CREATE TABLE kmers (
    id integer NOT NULL,
    hash integer,
    annotation_id integer
);


ALTER TABLE kmers OWNER TO gideon;

--
-- Name: kmers_id_seq; Type: SEQUENCE; Schema: public; Owner: gideon
--

CREATE SEQUENCE kmers_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE kmers_id_seq OWNER TO gideon;

--
-- Name: kmers_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gideon
--

ALTER SEQUENCE kmers_id_seq OWNED BY kmers.id;


--
-- Name: tags; Type: TABLE; Schema: public; Owner: gideon; Tablespace: 
--

CREATE TABLE tags (
    id integer NOT NULL,
    tag text
);


ALTER TABLE tags OWNER TO gideon;

--
-- Name: tags_id_seq; Type: SEQUENCE; Schema: public; Owner: gideon
--

CREATE SEQUENCE tags_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE tags_id_seq OWNER TO gideon;

--
-- Name: tags_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: gideon
--

ALTER SEQUENCE tags_id_seq OWNED BY tags.id;


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: gideon
--

ALTER TABLE ONLY annotations ALTER COLUMN id SET DEFAULT nextval('annotations_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: gideon
--

ALTER TABLE ONLY kmers ALTER COLUMN id SET DEFAULT nextval('kmers_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: gideon
--

ALTER TABLE ONLY tags ALTER COLUMN id SET DEFAULT nextval('tags_id_seq'::regclass);


--
-- Data for Name: annotations; Type: TABLE DATA; Schema: public; Owner: gideon
--

COPY annotations (id, tags) FROM stdin;
0	{0,1,2}
1	{0,1}
\.


--
-- Name: annotations_id_seq; Type: SEQUENCE SET; Schema: public; Owner: gideon
--

SELECT pg_catalog.setval('annotations_id_seq', 1, false);


--
-- Data for Name: kmers; Type: TABLE DATA; Schema: public; Owner: gideon
--

COPY kmers (id, hash, annotation_id) FROM stdin;
4	0	0
5	1	0
6	2	0
\.


--
-- Name: kmers_id_seq; Type: SEQUENCE SET; Schema: public; Owner: gideon
--

SELECT pg_catalog.setval('kmers_id_seq', 6, true);


--
-- Data for Name: tags; Type: TABLE DATA; Schema: public; Owner: gideon
--

COPY tags (id, tag) FROM stdin;
0	tag0
1	tag1
2	tag2
\.


--
-- Name: tags_id_seq; Type: SEQUENCE SET; Schema: public; Owner: gideon
--

SELECT pg_catalog.setval('tags_id_seq', 1, false);


--
-- Name: annotations_pkey; Type: CONSTRAINT; Schema: public; Owner: gideon; Tablespace: 
--

ALTER TABLE ONLY annotations
    ADD CONSTRAINT annotations_pkey PRIMARY KEY (id);


--
-- Name: kmers_pkey; Type: CONSTRAINT; Schema: public; Owner: gideon; Tablespace: 
--

ALTER TABLE ONLY kmers
    ADD CONSTRAINT kmers_pkey PRIMARY KEY (id);


--
-- Name: tags_pkey; Type: CONSTRAINT; Schema: public; Owner: gideon; Tablespace: 
--

ALTER TABLE ONLY tags
    ADD CONSTRAINT tags_pkey PRIMARY KEY (id);


--
-- Name: tags_tag_key; Type: CONSTRAINT; Schema: public; Owner: gideon; Tablespace: 
--

ALTER TABLE ONLY tags
    ADD CONSTRAINT tags_tag_key UNIQUE (tag);


--
-- Name: kmers_annotation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: gideon
--

ALTER TABLE ONLY kmers
    ADD CONSTRAINT kmers_annotation_id_fkey FOREIGN KEY (annotation_id) REFERENCES annotations(id);


--
-- Name: public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- PostgreSQL database dump complete
--

