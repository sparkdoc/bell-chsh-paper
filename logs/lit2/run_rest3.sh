#!/bin/bash
fetch() {
  for attempt in 1 2 3; do
    code=$(curl -sL --max-time 90 -w "%{http_code}" "https://export.arxiv.org/api/query?search_query=${2}&sortBy=submittedDate&sortOrder=descending&max_results=20" -o "$1.xml")
    echo "$(date -u +%H:%M:%S) $1 attempt$attempt HTTP $code size $(wc -c < $1.xml)" >> status3.log
    if [ "$code" = "200" ]; then return 0; fi
    sleep 90
  done
  echo "$(date -u +%H:%M:%S) $1 FAILED after 3 attempts" >> status3.log
}
fetch q11 'au:%22Alai%2C+Aaron%22+AND+cat:quant-ph'
sleep 60
fetch q12 'au:%22Pal%2C+Prosanta%22'
echo DONE >> status3.log
