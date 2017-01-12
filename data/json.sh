export PATH="/home/atp/kazama/data1/app/bin:$PATH"

cat dummy.json | jq-linux64 '.'  > pax_metadata.json
rm dummy.json
