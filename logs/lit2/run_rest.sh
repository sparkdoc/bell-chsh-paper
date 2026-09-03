#!/bin/bash
# Retry remaining arXiv API queries with long spacing (API is rate-limiting this IP)
fetch() { # $1=name $2=query
  for attempt in 1 2 3; do
    code=$(curl -sL --max-time 90 -w "%{http_code}" "https://export.arxiv.org/api/query?search_query=${2}&sortBy=submittedDate&sortOrder=descending&max_results=50" -o "$1.xml")
    echo "$(date -u +%H:%M:%S) $1 attempt$attempt HTTP $code size $(wc -c < $1.xml)" >> status.log
    if [ "$code" = "200" ]; then return 0; fi
    sleep 90
  done
  echo "$(date -u +%H:%M:%S) $1 FAILED after 3 attempts" >> status.log
}
sleep 60
fetch q3 'all:%22fine-tuning%22+AND+all:%22Bell%22'
sleep 90
fetch q4 'all:%22statistical+independence+loophole%22'
sleep 90
fetch q5 'all:%22measurement+independence%22+AND+all:%22CHSH%22'
sleep 90
fetch q6 'all:%22no-signalling%22+AND+all:%22hidden+variable%22+AND+all:%22fine-tuning%22'
echo DONE >> status.log
