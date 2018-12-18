--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
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


--
-- Name: citation_kind_enum; Type: TYPE; Schema: public; Owner: postgres
--

CREATE TYPE public.citation_kind_enum AS ENUM (
    'EXTERNAL',
    'INTERNAL'
);


ALTER TYPE public.citation_kind_enum OWNER TO postgres;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: affiliations; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.affiliations (
    a_id integer NOT NULL,
    name character varying NOT NULL,
    address character varying NOT NULL,
    email_domain character varying
);


ALTER TABLE public.affiliations OWNER TO postgres;

--
-- Name: affiliations_a_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.affiliations_a_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.affiliations_a_id_seq OWNER TO postgres;

--
-- Name: affiliations_a_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.affiliations_a_id_seq OWNED BY public.affiliations.a_id;


--
-- Name: alembic_version; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.alembic_version (
    version_num character varying(32) NOT NULL
);


ALTER TABLE public.alembic_version OWNER TO postgres;

--
-- Name: atoms; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.atoms (
    atom_id integer NOT NULL,
    struct_id integer NOT NULL,
    number integer NOT NULL,
    x double precision NOT NULL,
    y double precision NOT NULL,
    z double precision NOT NULL,
    charge double precision,
    magmom double precision,
    rmt double precision
);


ALTER TABLE public.atoms OWNER TO postgres;

--
-- Name: atoms_atom_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.atoms_atom_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.atoms_atom_id_seq OWNER TO postgres;

--
-- Name: atoms_atom_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.atoms_atom_id_seq OWNED BY public.atoms.atom_id;


--
-- Name: basis_sets; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.basis_sets (
    calc_id integer NOT NULL,
    type character varying NOT NULL,
    rgkmax double precision,
    lmaxapw double precision,
    lmaxmat double precision,
    lmaxvr double precision,
    gmaxvr double precision,
    repr bytea
);


ALTER TABLE public.basis_sets OWNER TO postgres;

--
-- Name: calcsets; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.calcsets (
    parent_calc_id integer NOT NULL,
    children_calc_id integer NOT NULL
);


ALTER TABLE public.calcsets OWNER TO postgres;

--
-- Name: calculations; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.calculations (
    checksum character varying,
    siblings_count integer,
    pottype_id integer,
    origin_id integer,
    nested_depth integer,
    frozen boolean,
    calc_id integer NOT NULL
);


ALTER TABLE public.calculations OWNER TO postgres;

--
-- Name: calculations_calc_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.calculations_calc_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.calculations_calc_id_seq OWNER TO postgres;

--
-- Name: calculations_calc_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.calculations_calc_id_seq OWNED BY public.calculations.calc_id;


--
-- Name: charges; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.charges (
    calc_id integer NOT NULL,
    core double precision NOT NULL,
    leakage double precision NOT NULL,
    valence double precision NOT NULL,
    interstitial double precision NOT NULL,
    muffintins double precision NOT NULL,
    total double precision NOT NULL
);


ALTER TABLE public.charges OWNER TO postgres;

--
-- Name: citations; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.citations (
    citation_id integer NOT NULL,
    value character varying NOT NULL,
    kind public.citation_kind_enum
);


ALTER TABLE public.citations OWNER TO postgres;

--
-- Name: citations_citation_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.citations_citation_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.citations_citation_id_seq OWNER TO postgres;

--
-- Name: citations_citation_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.citations_citation_id_seq OWNED BY public.citations.citation_id;


--
-- Name: coauthorships; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.coauthorships (
    calc_id integer,
    user_id integer
);


ALTER TABLE public.coauthorships OWNER TO postgres;

--
-- Name: codefamilies; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.codefamilies (
    family_id integer NOT NULL,
    content character varying NOT NULL
);


ALTER TABLE public.codefamilies OWNER TO postgres;

--
-- Name: codefamilies_family_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.codefamilies_family_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.codefamilies_family_id_seq OWNER TO postgres;

--
-- Name: codefamilies_family_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.codefamilies_family_id_seq OWNED BY public.codefamilies.family_id;


--
-- Name: codeversions; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.codeversions (
    version_id integer NOT NULL,
    family_id integer,
    content character varying NOT NULL
);


ALTER TABLE public.codeversions OWNER TO postgres;

--
-- Name: codeversions_version_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.codeversions_version_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.codeversions_version_id_seq OWNER TO postgres;

--
-- Name: codeversions_version_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.codeversions_version_id_seq OWNED BY public.codeversions.version_id;


--
-- Name: doi_mapping; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.doi_mapping (
    calc_id integer NOT NULL,
    id_str character varying NOT NULL
);


ALTER TABLE public.doi_mapping OWNER TO postgres;

--
-- Name: doi_mapping_calc_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.doi_mapping_calc_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.doi_mapping_calc_id_seq OWNER TO postgres;

--
-- Name: doi_mapping_calc_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.doi_mapping_calc_id_seq OWNED BY public.doi_mapping.calc_id;


--
-- Name: eigenvalues; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.eigenvalues (
    eid integer NOT NULL,
    electrons_calc_id integer,
    phonons_calc_id integer,
    dos bytea,
    bands bytea,
    projected bytea,
    eigenvalues bytea
);


ALTER TABLE public.eigenvalues OWNER TO postgres;

--
-- Name: eigenvalues_eid_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.eigenvalues_eid_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.eigenvalues_eid_seq OWNER TO postgres;

--
-- Name: eigenvalues_eid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.eigenvalues_eid_seq OWNED BY public.eigenvalues.eid;


--
-- Name: electrons; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.electrons (
    calc_id integer NOT NULL,
    gap double precision,
    is_direct integer
);


ALTER TABLE public.electrons OWNER TO postgres;

--
-- Name: energies; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.energies (
    calc_id integer NOT NULL,
    convergence bytea,
    total double precision
);


ALTER TABLE public.energies OWNER TO postgres;

--
-- Name: forces; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.forces (
    calc_id integer NOT NULL,
    "values" bytea NOT NULL
);


ALTER TABLE public.forces OWNER TO postgres;

--
-- Name: grid; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.grid (
    calc_id integer NOT NULL,
    info bytea
);


ALTER TABLE public.grid OWNER TO postgres;

--
-- Name: lattices; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.lattices (
    struct_id integer NOT NULL,
    a double precision NOT NULL,
    b double precision NOT NULL,
    c double precision NOT NULL,
    alpha double precision NOT NULL,
    beta double precision NOT NULL,
    gamma double precision NOT NULL,
    a11 double precision NOT NULL,
    a12 double precision NOT NULL,
    a13 double precision NOT NULL,
    a21 double precision NOT NULL,
    a22 double precision NOT NULL,
    a23 double precision NOT NULL,
    a31 double precision NOT NULL,
    a32 double precision NOT NULL,
    a33 double precision NOT NULL
);


ALTER TABLE public.lattices OWNER TO postgres;

--
-- Name: login_tokens; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.login_tokens (
    user_id integer,
    token character varying NOT NULL,
    valid_until timestamp with time zone
);


ALTER TABLE public.login_tokens OWNER TO postgres;

--
-- Name: metadata; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.metadata (
    calc_id integer NOT NULL,
    version_id integer,
    location character varying,
    finished integer,
    raw_input text,
    modeling_time double precision,
    chemical_formula character varying,
    added timestamp with time zone,
    oadate timestamp with time zone,
    download_size bigint,
    filenames bytea
);


ALTER TABLE public.metadata OWNER TO postgres;

--
-- Name: metadata_citations; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.metadata_citations (
    calc_id integer,
    citation_id integer
);


ALTER TABLE public.metadata_citations OWNER TO postgres;

--
-- Name: ownerships; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.ownerships (
    calc_id integer,
    user_id integer
);


ALTER TABLE public.ownerships OWNER TO postgres;

--
-- Name: phonons; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.phonons (
    calc_id integer NOT NULL
);


ALTER TABLE public.phonons OWNER TO postgres;

--
-- Name: pottypes; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.pottypes (
    pottype_id integer NOT NULL,
    name character varying NOT NULL
);


ALTER TABLE public.pottypes OWNER TO postgres;

--
-- Name: pottypes_pottype_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.pottypes_pottype_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.pottypes_pottype_id_seq OWNER TO postgres;

--
-- Name: pottypes_pottype_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.pottypes_pottype_id_seq OWNED BY public.pottypes.pottype_id;


--
-- Name: pragma; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.pragma (
    content character varying NOT NULL
);


ALTER TABLE public.pragma OWNER TO postgres;

--
-- Name: recipintegs; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.recipintegs (
    calc_id integer NOT NULL,
    kgrid character varying,
    kshift double precision,
    smearing double precision,
    smeartype character varying
);


ALTER TABLE public.recipintegs OWNER TO postgres;

--
-- Name: sessions; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.sessions (
    token character varying NOT NULL,
    user_id integer,
    valid_until timestamp with time zone,
    last_access timestamp with time zone,
    permission integer
);


ALTER TABLE public.sessions OWNER TO postgres;

--
-- Name: shareships; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.shareships (
    calc_id integer,
    user_id integer
);


ALTER TABLE public.shareships OWNER TO postgres;

--
-- Name: spacegroups; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.spacegroups (
    calc_id integer NOT NULL,
    n integer NOT NULL
);


ALTER TABLE public.spacegroups OWNER TO postgres;

--
-- Name: struct_optimisation; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.struct_optimisation (
    calc_id integer NOT NULL,
    tresholds bytea,
    ncycles bytea
);


ALTER TABLE public.struct_optimisation OWNER TO postgres;

--
-- Name: struct_ratios; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.struct_ratios (
    calc_id integer NOT NULL,
    chemical_formula character varying NOT NULL,
    is_primitive boolean,
    formula_units integer NOT NULL,
    nelem integer NOT NULL,
    dimensions double precision
);


ALTER TABLE public.struct_ratios OWNER TO postgres;

--
-- Name: structures; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.structures (
    struct_id integer NOT NULL,
    calc_id integer NOT NULL,
    step integer NOT NULL,
    final boolean NOT NULL
);


ALTER TABLE public.structures OWNER TO postgres;

--
-- Name: structures_struct_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.structures_struct_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.structures_struct_id_seq OWNER TO postgres;

--
-- Name: structures_struct_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.structures_struct_id_seq OWNED BY public.structures.struct_id;


--
-- Name: tags; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.tags (
    calc_id integer NOT NULL,
    tid integer NOT NULL
);


ALTER TABLE public.tags OWNER TO postgres;

--
-- Name: topics; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.topics (
    tid integer NOT NULL,
    cid integer NOT NULL,
    topic character varying
);


ALTER TABLE public.topics OWNER TO postgres;

--
-- Name: topics_tid_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.topics_tid_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.topics_tid_seq OWNER TO postgres;

--
-- Name: topics_tid_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.topics_tid_seq OWNED BY public.topics.tid;


--
-- Name: uploads; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.uploads (
    upload_id integer NOT NULL,
    upload_name character varying NOT NULL,
    user_id integer,
    is_processed boolean NOT NULL,
    created timestamp with time zone,
    is_all_uploaded boolean,
    target_path character varying,
    skip_extraction boolean
);


ALTER TABLE public.uploads OWNER TO postgres;

--
-- Name: uploads_upload_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.uploads_upload_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.uploads_upload_id_seq OWNER TO postgres;

--
-- Name: uploads_upload_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.uploads_upload_id_seq OWNED BY public.uploads.upload_id;


--
-- Name: user_metadata; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.user_metadata (
    calc_id integer NOT NULL,
    permission integer,
    label character varying
);


ALTER TABLE public.user_metadata OWNER TO postgres;

--
-- Name: users; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE public.users (
    user_id integer NOT NULL,
    firstname character varying NOT NULL,
    lastname character varying NOT NULL,
    username character varying,
    email character varying,
    affiliation integer,
    password character varying,
    created timestamp with time zone
);


ALTER TABLE public.users OWNER TO postgres;

--
-- Name: users_user_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.users_user_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.users_user_id_seq OWNER TO postgres;

--
-- Name: users_user_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.users_user_id_seq OWNED BY public.users.user_id;


--
-- Name: a_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.affiliations ALTER COLUMN a_id SET DEFAULT nextval('public.affiliations_a_id_seq'::regclass);


--
-- Name: atom_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.atoms ALTER COLUMN atom_id SET DEFAULT nextval('public.atoms_atom_id_seq'::regclass);


--
-- Name: calc_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.calculations ALTER COLUMN calc_id SET DEFAULT nextval('public.calculations_calc_id_seq'::regclass);


--
-- Name: citation_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.citations ALTER COLUMN citation_id SET DEFAULT nextval('public.citations_citation_id_seq'::regclass);


--
-- Name: family_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.codefamilies ALTER COLUMN family_id SET DEFAULT nextval('public.codefamilies_family_id_seq'::regclass);


--
-- Name: version_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.codeversions ALTER COLUMN version_id SET DEFAULT nextval('public.codeversions_version_id_seq'::regclass);


--
-- Name: calc_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.doi_mapping ALTER COLUMN calc_id SET DEFAULT nextval('public.doi_mapping_calc_id_seq'::regclass);


--
-- Name: eid; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.eigenvalues ALTER COLUMN eid SET DEFAULT nextval('public.eigenvalues_eid_seq'::regclass);


--
-- Name: pottype_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pottypes ALTER COLUMN pottype_id SET DEFAULT nextval('public.pottypes_pottype_id_seq'::regclass);


--
-- Name: struct_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.structures ALTER COLUMN struct_id SET DEFAULT nextval('public.structures_struct_id_seq'::regclass);


--
-- Name: tid; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.topics ALTER COLUMN tid SET DEFAULT nextval('public.topics_tid_seq'::regclass);


--
-- Name: upload_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.uploads ALTER COLUMN upload_id SET DEFAULT nextval('public.uploads_upload_id_seq'::regclass);


--
-- Name: user_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.users ALTER COLUMN user_id SET DEFAULT nextval('public.users_user_id_seq'::regclass);


--
-- Data for Name: affiliations; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: affiliations_a_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.affiliations_a_id_seq', 1, false);


--
-- Data for Name: alembic_version; Type: TABLE DATA; Schema: public; Owner: postgres
--

INSERT INTO public.alembic_version VALUES ('12ed242720e1');


--
-- Data for Name: atoms; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: atoms_atom_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.atoms_atom_id_seq', 1, false);


--
-- Data for Name: basis_sets; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: calcsets; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: calculations; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: calculations_calc_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.calculations_calc_id_seq', 1, false);


--
-- Data for Name: charges; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: citations; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: citations_citation_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.citations_citation_id_seq', 1, false);


--
-- Data for Name: coauthorships; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: codefamilies; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: codefamilies_family_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.codefamilies_family_id_seq', 1, false);


--
-- Data for Name: codeversions; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: codeversions_version_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.codeversions_version_id_seq', 1, false);


--
-- Data for Name: doi_mapping; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: doi_mapping_calc_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.doi_mapping_calc_id_seq', 1, false);


--
-- Data for Name: eigenvalues; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: eigenvalues_eid_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.eigenvalues_eid_seq', 1, false);


--
-- Data for Name: electrons; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: energies; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: forces; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: grid; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: lattices; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: login_tokens; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: metadata; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: metadata_citations; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: ownerships; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: phonons; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: pottypes; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: pottypes_pottype_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.pottypes_pottype_id_seq', 1, false);


--
-- Data for Name: pragma; Type: TABLE DATA; Schema: public; Owner: postgres
--

INSERT INTO public.pragma VALUES ('4.59');


--
-- Data for Name: recipintegs; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: sessions; Type: TABLE DATA; Schema: public; Owner: postgres
--

INSERT INTO public.sessions VALUES ('leonard.hofstadter@nomad-fairdi.tests.de', 2, '2100-12-17 09:00:00+00', NULL, NULL);
INSERT INTO public.sessions VALUES ('sheldon.cooper@nomad-fairdi.tests.de', 1, '2100-12-17 09:00:00+00', NULL, NULL);


--
-- Data for Name: shareships; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: spacegroups; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: struct_optimisation; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: struct_ratios; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: structures; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: structures_struct_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.structures_struct_id_seq', 1, false);


--
-- Data for Name: tags; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: topics; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: topics_tid_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.topics_tid_seq', 1, false);


--
-- Data for Name: uploads; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Name: uploads_upload_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.uploads_upload_id_seq', 1, false);


--
-- Data for Name: user_metadata; Type: TABLE DATA; Schema: public; Owner: postgres
--



--
-- Data for Name: users; Type: TABLE DATA; Schema: public; Owner: postgres
--

INSERT INTO public.users VALUES (1, 'Sheldon', 'Cooper', 'sheldon.cooper', 'sheldon.cooper@nomad-fairdi.tests.de', NULL, '$2y$12$jths1LQPsLofuBQ3evVIluhQeQ/BZfbdTSZHFcPGdcNmHz2WvDj.y', NULL);
INSERT INTO public.users VALUES (2, 'Leonard', 'Hofstadter', 'leonard.hofstadter', 'leonard.hofstadter@nomad-fairdi.tests.de', NULL, '$2y$12$jths1LQPsLofuBQ3evVIluhQeQ/BZfbdTSZHFcPGdcNmHz2WvDj.y', NULL);


--
-- Name: users_user_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.users_user_id_seq', 1, true);


--
-- Name: affiliations_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.affiliations
    ADD CONSTRAINT affiliations_pkey PRIMARY KEY (a_id);


--
-- Name: atoms_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.atoms
    ADD CONSTRAINT atoms_pkey PRIMARY KEY (atom_id);


--
-- Name: basis_sets_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.basis_sets
    ADD CONSTRAINT basis_sets_pkey PRIMARY KEY (calc_id);


--
-- Name: calculations_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.calculations
    ADD CONSTRAINT calculations_pkey PRIMARY KEY (calc_id);


--
-- Name: charges_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.charges
    ADD CONSTRAINT charges_pkey PRIMARY KEY (calc_id);


--
-- Name: citations_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.citations
    ADD CONSTRAINT citations_pkey PRIMARY KEY (citation_id);


--
-- Name: citations_value_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.citations
    ADD CONSTRAINT citations_value_key UNIQUE (value);


--
-- Name: codefamilies_content_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.codefamilies
    ADD CONSTRAINT codefamilies_content_key UNIQUE (content);


--
-- Name: codefamilies_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.codefamilies
    ADD CONSTRAINT codefamilies_pkey PRIMARY KEY (family_id);


--
-- Name: codeversions_content_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.codeversions
    ADD CONSTRAINT codeversions_content_key UNIQUE (content);


--
-- Name: codeversions_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.codeversions
    ADD CONSTRAINT codeversions_pkey PRIMARY KEY (version_id);


--
-- Name: doi_mapping_id_str_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.doi_mapping
    ADD CONSTRAINT doi_mapping_id_str_key UNIQUE (id_str);


--
-- Name: doi_mapping_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.doi_mapping
    ADD CONSTRAINT doi_mapping_pkey PRIMARY KEY (calc_id);


--
-- Name: eigenvalues_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.eigenvalues
    ADD CONSTRAINT eigenvalues_pkey PRIMARY KEY (eid);


--
-- Name: electrons_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.electrons
    ADD CONSTRAINT electrons_pkey PRIMARY KEY (calc_id);


--
-- Name: energies_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.energies
    ADD CONSTRAINT energies_pkey PRIMARY KEY (calc_id);


--
-- Name: forces_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.forces
    ADD CONSTRAINT forces_pkey PRIMARY KEY (calc_id);


--
-- Name: grid_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.grid
    ADD CONSTRAINT grid_pkey PRIMARY KEY (calc_id);


--
-- Name: lattices_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.lattices
    ADD CONSTRAINT lattices_pkey PRIMARY KEY (struct_id);


--
-- Name: login_tokens_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.login_tokens
    ADD CONSTRAINT login_tokens_pkey PRIMARY KEY (token);


--
-- Name: metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.metadata
    ADD CONSTRAINT metadata_pkey PRIMARY KEY (calc_id);


--
-- Name: phonons_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.phonons
    ADD CONSTRAINT phonons_pkey PRIMARY KEY (calc_id);


--
-- Name: pottypes_name_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.pottypes
    ADD CONSTRAINT pottypes_name_key UNIQUE (name);


--
-- Name: pottypes_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.pottypes
    ADD CONSTRAINT pottypes_pkey PRIMARY KEY (pottype_id);


--
-- Name: pragma_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.pragma
    ADD CONSTRAINT pragma_pkey PRIMARY KEY (content);


--
-- Name: recipintegs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.recipintegs
    ADD CONSTRAINT recipintegs_pkey PRIMARY KEY (calc_id);


--
-- Name: sessions_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.sessions
    ADD CONSTRAINT sessions_pkey PRIMARY KEY (token);


--
-- Name: spacegroups_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.spacegroups
    ADD CONSTRAINT spacegroups_pkey PRIMARY KEY (calc_id);


--
-- Name: struct_optimisation_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.struct_optimisation
    ADD CONSTRAINT struct_optimisation_pkey PRIMARY KEY (calc_id);


--
-- Name: struct_ratios_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.struct_ratios
    ADD CONSTRAINT struct_ratios_pkey PRIMARY KEY (calc_id);


--
-- Name: structures_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.structures
    ADD CONSTRAINT structures_pkey PRIMARY KEY (struct_id);


--
-- Name: topics_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.topics
    ADD CONSTRAINT topics_pkey PRIMARY KEY (tid);


--
-- Name: u_children_parent_calc_id; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.calcsets
    ADD CONSTRAINT u_children_parent_calc_id UNIQUE (children_calc_id, parent_calc_id);


--
-- Name: u_coauthorships_calc_id_user; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.coauthorships
    ADD CONSTRAINT u_coauthorships_calc_id_user UNIQUE (calc_id, user_id);


--
-- Name: u_coauthorships_user_calc_id; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.coauthorships
    ADD CONSTRAINT u_coauthorships_user_calc_id UNIQUE (user_id, calc_id);


--
-- Name: u_metadata_citations_calc_citation; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.metadata_citations
    ADD CONSTRAINT u_metadata_citations_calc_citation UNIQUE (calc_id, citation_id);


--
-- Name: u_metadata_citations_citation_calc; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.metadata_citations
    ADD CONSTRAINT u_metadata_citations_citation_calc UNIQUE (citation_id, calc_id);


--
-- Name: u_ownerships_calc_id_user; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.ownerships
    ADD CONSTRAINT u_ownerships_calc_id_user UNIQUE (calc_id, user_id);


--
-- Name: u_ownerships_user_calc_id; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.ownerships
    ADD CONSTRAINT u_ownerships_user_calc_id UNIQUE (user_id, calc_id);


--
-- Name: u_parent_children_calc_id; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.calcsets
    ADD CONSTRAINT u_parent_children_calc_id PRIMARY KEY (parent_calc_id, children_calc_id);


--
-- Name: u_shareships_calc_id_user; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.shareships
    ADD CONSTRAINT u_shareships_calc_id_user UNIQUE (calc_id, user_id);


--
-- Name: u_shareships_user_calc_id; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.shareships
    ADD CONSTRAINT u_shareships_user_calc_id UNIQUE (user_id, calc_id);


--
-- Name: u_tags_calc_id_tid; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT u_tags_calc_id_tid UNIQUE (calc_id, tid);


--
-- Name: u_tags_tid_calc_id; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT u_tags_tid_calc_id UNIQUE (tid, calc_id);


--
-- Name: uploads_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.uploads
    ADD CONSTRAINT uploads_pkey PRIMARY KEY (upload_id);


--
-- Name: user_metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.user_metadata
    ADD CONSTRAINT user_metadata_pkey PRIMARY KEY (calc_id);


--
-- Name: users_email_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_email_key UNIQUE (email);


--
-- Name: users_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_pkey PRIMARY KEY (user_id);


--
-- Name: users_username_key; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_username_key UNIQUE (username);


--
-- Name: atoms_struct_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.atoms
    ADD CONSTRAINT atoms_struct_id_fkey FOREIGN KEY (struct_id) REFERENCES public.structures(struct_id);


--
-- Name: basis_sets_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.basis_sets
    ADD CONSTRAINT basis_sets_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: calcsets_children_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.calcsets
    ADD CONSTRAINT calcsets_children_calc_id_fkey FOREIGN KEY (children_calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: calcsets_parent_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.calcsets
    ADD CONSTRAINT calcsets_parent_calc_id_fkey FOREIGN KEY (parent_calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: calculations_origin_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.calculations
    ADD CONSTRAINT calculations_origin_id_fkey FOREIGN KEY (origin_id) REFERENCES public.uploads(upload_id);


--
-- Name: calculations_pottype_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.calculations
    ADD CONSTRAINT calculations_pottype_id_fkey FOREIGN KEY (pottype_id) REFERENCES public.pottypes(pottype_id);


--
-- Name: charges_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.charges
    ADD CONSTRAINT charges_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: coauthorships_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.coauthorships
    ADD CONSTRAINT coauthorships_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: coauthorships_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.coauthorships
    ADD CONSTRAINT coauthorships_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(user_id);


--
-- Name: codeversions_family_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.codeversions
    ADD CONSTRAINT codeversions_family_id_fkey FOREIGN KEY (family_id) REFERENCES public.codefamilies(family_id);


--
-- Name: eigenvalues_electrons_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.eigenvalues
    ADD CONSTRAINT eigenvalues_electrons_calc_id_fkey FOREIGN KEY (electrons_calc_id) REFERENCES public.electrons(calc_id);


--
-- Name: eigenvalues_phonons_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.eigenvalues
    ADD CONSTRAINT eigenvalues_phonons_calc_id_fkey FOREIGN KEY (phonons_calc_id) REFERENCES public.phonons(calc_id);


--
-- Name: electrons_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.electrons
    ADD CONSTRAINT electrons_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: energies_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.energies
    ADD CONSTRAINT energies_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: forces_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.forces
    ADD CONSTRAINT forces_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: grid_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.grid
    ADD CONSTRAINT grid_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: lattices_struct_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.lattices
    ADD CONSTRAINT lattices_struct_id_fkey FOREIGN KEY (struct_id) REFERENCES public.structures(struct_id);


--
-- Name: login_tokens_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.login_tokens
    ADD CONSTRAINT login_tokens_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(user_id);


--
-- Name: metadata_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.metadata
    ADD CONSTRAINT metadata_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: metadata_citations_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.metadata_citations
    ADD CONSTRAINT metadata_citations_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: metadata_citations_citation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.metadata_citations
    ADD CONSTRAINT metadata_citations_citation_id_fkey FOREIGN KEY (citation_id) REFERENCES public.citations(citation_id);


--
-- Name: metadata_version_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.metadata
    ADD CONSTRAINT metadata_version_id_fkey FOREIGN KEY (version_id) REFERENCES public.codeversions(version_id);


--
-- Name: ownerships_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.ownerships
    ADD CONSTRAINT ownerships_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: ownerships_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.ownerships
    ADD CONSTRAINT ownerships_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(user_id);


--
-- Name: phonons_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.phonons
    ADD CONSTRAINT phonons_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: recipintegs_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.recipintegs
    ADD CONSTRAINT recipintegs_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: sessions_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.sessions
    ADD CONSTRAINT sessions_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(user_id);


--
-- Name: shareships_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.shareships
    ADD CONSTRAINT shareships_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: shareships_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.shareships
    ADD CONSTRAINT shareships_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(user_id);


--
-- Name: spacegroups_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.spacegroups
    ADD CONSTRAINT spacegroups_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: struct_optimisation_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.struct_optimisation
    ADD CONSTRAINT struct_optimisation_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: struct_ratios_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.struct_ratios
    ADD CONSTRAINT struct_ratios_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: structures_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.structures
    ADD CONSTRAINT structures_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: tags_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT tags_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: tags_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT tags_tid_fkey FOREIGN KEY (tid) REFERENCES public.topics(tid);


--
-- Name: uploads_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.uploads
    ADD CONSTRAINT uploads_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(user_id);


--
-- Name: user_metadata_calc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.user_metadata
    ADD CONSTRAINT user_metadata_calc_id_fkey FOREIGN KEY (calc_id) REFERENCES public.calculations(calc_id);


--
-- Name: users_affiliation_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_affiliation_fkey FOREIGN KEY (affiliation) REFERENCES public.affiliations(a_id);


--
-- Name: SCHEMA public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- PostgreSQL database dump complete
--

