s/^\(DROP INDEX \)/\1IF EXISTS /;
s/^\(DROP TABLE \)/\1IF EXISTS /;
s/ NONE / BLOB /;
s/'NULL'/NULL/;
s///;
