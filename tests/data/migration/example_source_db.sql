SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET client_min_messages = warning;

TRUNCATE TABLE public.users CASCADE;
INSERT INTO public.users VALUES (3, 'one', 'one', 'one', 'one', NULL, '$2y$12$jths1LQPsLofuBQ3evVIluhQeQ/BZfbdTSZHFcPGdcNmHz2WvDj.y', NULL);
INSERT INTO public.users VALUES (4, 'two', 'two', 'two', 'two', NULL, '$2y$12$jths1LQPsLofuBQ3evVIluhQeQ/BZfbdTSZHFcPGdcNmHz2WvDj.y', NULL);
INSERT INTO public.calculations VALUES (NULL, NULL, NULL, NULL, 0, false, 1, NULL);
INSERT INTO public.calculations VALUES (NULL, NULL, NULL, NULL, 0, false, 2, NULL);
INSERT INTO public.codefamilies VALUES (1, 'VASP');
INSERT INTO public.codeversions VALUES (1, 1, '4.6.35');
-- topcis
INSERT INTO public.topics VALUES (1, 90, 'tetragonal');
INSERT INTO public.topics VALUES (2, 220, 'VASP');
INSERT INTO public.topics VALUES (3, 50, 'bulk');
INSERT INTO public.topics VALUES (4, 75, 'GGA');
INSERT INTO public.topics VALUES (5, 80, 'plane waves');
INSERT INTO public.topics VALUES (6, 10, 'Br');
INSERT INTO public.topics VALUES (7, 10, 'K');
INSERT INTO public.topics VALUES (8, 10, 'Si');
-- mapping topics to calcs via tags
INSERT INTO public.tags VALUES(1, 1);
INSERT INTO public.tags VALUES(2, 1);
INSERT INTO public.tags VALUES(1, 2);
INSERT INTO public.tags VALUES(2, 2);
INSERT INTO public.tags VALUES(1, 3);
INSERT INTO public.tags VALUES(2, 3);
INSERT INTO public.tags VALUES(1, 4);
INSERT INTO public.tags VALUES(2, 4);
INSERT INTO public.tags VALUES(1, 5);
INSERT INTO public.tags VALUES(2, 5);
INSERT INTO public.tags VALUES(1, 6);
INSERT INTO public.tags VALUES(2, 6);
INSERT INTO public.tags VALUES(1, 7);
INSERT INTO public.tags VALUES(2, 7);
INSERT INTO public.tags VALUES(1, 8);
INSERT INTO public.tags VALUES(2, 8);

INSERT INTO public.metadata VALUES (1, NULL, NULL, NULL, NULL, 'BrKSi2', '2019-01-01 12:00:00', NULL, decode('["$EXTRACTED/upload/1/template.json"]', 'escape'), 1, NULL);
INSERT INTO public.metadata VALUES (1, NULL, NULL, NULL, NULL, 'BrKSi2', '2015-01-01 13:00:00', NULL, decode('["$EXTRACTED/upload/2/template.json"]', 'escape'), 2, NULL);
INSERT INTO public.spacegroups VALUES (1, 123);
INSERT INTO public.spacegroups VALUES (2, 123);
INSERT INTO public.user_metadata VALUES (1, 0, 'label1');
INSERT INTO public.user_metadata VALUES (2, 1, 'label2');
INSERT INTO public.ownerships VALUES (1, 3);
INSERT INTO public.ownerships VALUES (2, 4);
INSERT INTO public.coauthorships VALUES (1, 4);
INSERT INTO public.shareships VALUES (2, 3);

-- example dataset
INSERT INTO public.calculations VALUES (NULL, NULL, NULL, NULL, 1, false, 3, NULL);
INSERT INTO public.calcsets VALUES (3, 1);
INSERT INTO public.calcsets VALUES (3, 2);
INSERT INTO public.citations VALUES(1, 'internal_ref', 'INTERNAL');
INSERT INTO public.metadata_citations VALUES (3, 1);
INSERT INTO public.metadata VALUES (NULL, NULL, NULL, NULL, NULL, 'test_dataset', '2019-01-01 12:00:00', NULL, NULL, 3, NULL);

-- example ref
INSERT INTO public.citations VALUES(2, 'external_ref', 'EXTERNAL');
INSERT INTO public.metadata_citations VALUES (1, 2);