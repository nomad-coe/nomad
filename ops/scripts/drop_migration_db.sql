UPDATE pg_database SET datallowconn = 'false' WHERE datname = 'fairdi_nomad_migration';
SELECT pg_terminate_backend(pid)
FROM pg_stat_activity
WHERE datname = 'fairdi_nomad_migration';
DROP DATABASE fairdi_nomad_migration;
