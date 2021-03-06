#!/bin/sh

case $# in
0)
	exit 1
	;;
esac

file="$1"

# extract master
sed '
/^DROP INDEX IF EXISTS IDX_MAIN_/d;
/^DROP TABLE IF EXISTS MAIN;/d;
/^DROP TABLE IF EXISTS MAIN_COORDINATE;/d;
/^DROP TABLE IF EXISTS FLAG_CMD;/d;
/^DROP TABLE IF EXISTS HISTORY_COMMAND;/d;
/^DROP TABLE IF EXISTS HISTORY_PARAM;/d;
/^DROP TABLE IF EXISTS HISTORY;/d;
/^CREATE TABLE MAIN$/,/^);/d;
/^CREATE TABLE MAIN_COORDINATE$/,/^);/d;
/^CREATE TABLE FLAG_CMD$/,/^);/d;
/^CREATE TABLE HISTORY_COMMAND$/,/^);/d;
/^CREATE TABLE HISTORY_PARAM$/,/^);/d;
/^CREATE TABLE HISTORY$/,/^);/d;
/^CREATE INDEX IDX_MAIN_/d;
' "$file" > MSM.ddl

cat << EOF > MST.ddl
-- ATTACH DATABASE 'ms.mdb' AS msm;

EOF
sed -n '
/^\/\*/p;
/^DROP INDEX IF EXISTS IDX_MAIN_/p;
/^DROP TABLE IF EXISTS MAIN;/p;
/^DROP TABLE IF EXISTS MAIN_COORDINATE;/p;
/^DROP TABLE IF EXISTS FLAG_CMD;/p;
/^DROP TABLE IF EXISTS HISTORY_COMMAND;/p;
/^DROP TABLE IF EXISTS HISTORY_PARAM;/p;
/^DROP TABLE IF EXISTS HISTORY;/p;
/^DROP TABLE IF EXISTS COLUMN_KEYWORD;/p;
/^DROP TABLE IF EXISTS TABLE_KEYWORD;/p;
/^CREATE TABLE MAIN$/,/^);/p;
/^CREATE TABLE MAIN_COORDINATE$/,/^);/p;
/^CREATE TABLE FLAG_CMD$/,/^);/p;
/^CREATE TABLE HISTORY_COMMAND$/,/^);/p;
/^CREATE TABLE HISTORY_PARAM$/,/^);/p;
/^CREATE TABLE HISTORY$/,/^);/p;
/^CREATE TABLE COLUMN_KEYWORD$/,/^);/p;
/^CREATE TABLE TABLE_KEYWORD$/,/^);/p;
/^CREATE INDEX IDX_MAIN_/p;
' "$file" >> MST.ddl

