#!/bin/bash
fetch() { # $1=name $2=query $3=sortfield
  for attempt in 1 2 3; do
    code=$(curl -sL --max-time 90 -w "%{http_code}" "https://export.arxiv.org/api/query?search_query=${2}&sortBy=${3}&sortOrder=descending&max_results=50" -o "$1.xml")
    echo "$(date -u +%H:%M:%S) $1 attempt$attempt HTTP $code size $(wc -c < $1.xml)" >> status2.log
    if [ "$code" = "200" ]; then return 0; fi
    sleep 90
  done
  echo "$(date -u +%H:%M:%S) $1 FAILED after 3 attempts" >> status2.log
}
fetch q3 'all:%22fine-tuning%22+AND+all:%22Bell%22' submittedDate
sleep 90
fetch q7 'all:%22measurement+dependence%22' lastUpdatedDate
sleep 90
fetch q8 'all:%22superdeterminism%22' lastUpdatedDate
sleep 90
fetch q9 'jr:%22Foundations+of+Physics%22' submittedDate
sleep 90
fetch q10 'jr:%22Studies+in+History+and+Philosophy+of+Science%22' submittedDate
echo DONE >> status2.log
