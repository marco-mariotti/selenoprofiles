#! /usr/bin/gawk -f 
# normal output: (separated by space); names of target and query are normally just the first word of it

#target_name target_start target_end strand query_name query_start query_end target_aligned_seq query_aligned_seq evalue bit_score

#all positions are 1-based. strand is derived from "Frame " field and refers to target. target_start printed is always < target_end ; for query, to allow telling the strand in case of tblastx, query_start > query_end when query strand is "-"
 
#if REF_FILE is defined, filename start_line and end_line are printed in each line. if TLENGTH is defined, the target length is printed after the bits score. If QLENGTH is defined, the query length is printed after it; if FULL_QUERY or FULL_TARGET is ==1, the name displayed is the full name and not just the first word. if FULL is ==1, it is equivalent to set FULL_QUERY and FULL_TARGET to 1. NB: the full names are printed replacing white spaces with the character %, to maintain the columns as it FULL is not active.
BEGIN {start_pos=-1
q_start_pos=-1
line_start=-1
COMMENT=""
ALI=1
SCORE=1
if (FULL==1){ FULL_QUERY=1; FULL_TARGET=1}
}     

{
 

  ##3
  if ( /\(.*letters/)   {

	split($0, LINE_SPLIT, "letters")
	
	while (sub(/[\(\)QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm=\:\ ]/, "", LINE_SPLIT[1])){}

	query_length=LINE_SPLIT[1]

  in_query_line=0

  }

  ##2

  if (in_query_line && FULL_QUERY){
  query=query " " $0
 }

  ###1 

  if ( /^Query=/)   {
	  if (FULL_QUERY){ query=substr($0, 8) } #taking everything after "Query= "
	  else {query = $2}

    if (line_start==-1)
    {line_start = NR }
    in_query_line=1

  }



if ( /^Sbjct:/ )    
{  
	if (ALI){
		tline=$0
		sub(/Sbjct:/, "")
		while (sub(/[0-9 ]/, "")) {}
		subj_seq = subj_seq $0
		$0=tline 
	}

	while  ( sub(/[A-Z]|[a-z]|\*|\-|:|\./, " ") )	#    sub(/[:alpha:]|[:upper:]|:/, " ")   )
	{ 	 }
	if ( start_pos==-1 ){	start_pos=$1	}
	end_pos=$2

}
if ( /^Query:/ )    #:[:blank:]*[:digit:]+[:blank:]*[:alpha:]/)
{  
	if (ALI){
		tline=$0
		sub(/Query:/, "")
		while (sub(/[0-9 ]/, "")) {}
		query_seq = query_seq $0
		$0=tline 
	}

	while  ( sub(/[A-Z]|[a-z]|\*|\-|:|\./, " ") )	#    sub(/[:alpha:]|[:upper:]|:/, " ")   )
	{ 	 }
	if ( q_start_pos==-1 ){	q_start_pos=$1	}
	q_end_pos=$2
}

if (/^>/ || /^\#\#/ || /Score =/ || /^Reference:/ || /^Matrix:/)
{	
	if (q_start_pos != -1)	#already parsed one hit?
	{ stop_parsing=1}



}
if (stop_parsing==1){ #look down as well
	if (ALI){
		COMMENT=COMMENT " " subj_seq " " query_seq 
		
	}
	if (SCORE)
	{	COMMENT=COMMENT " " evalue " " bits 
	}

	if (TLENGTH)
	{       COMMENT=COMMENT " " target_length
	}

	if (QLENGTH)
	{	COMMENT=COMMENT " " query_length 
	}


	if (REF_FILE)
	{	COMMENT=COMMENT " " ARGV[1] " " line_start " " (NR-1) 

	}

	strand="+"
	if (start_pos>end_pos){ temp=start_pos; start_pos=end_pos; end_pos=temp; strand="-"} #exchanging start with stop. using n_hits as a temp var

	if (! query){		query="UNKNOWN_QUERY" }
	
	if (FULL_QUERY){	gsub(/ +/, "%", query) }
	if (FULL_TARGET){	gsub(/ +/, "%", subject) }	
	
	print subject " " start_pos " " end_pos " " strand " "  query " " q_start_pos " " q_end_pos COMMENT
	start_pos=-1
	end_pos=-1
	q_start_pos=-1
	q_end_pos=-1
	line_start=-1
	stop_parsing=0
	COMMENT=""
	subj_seq=""
	query_seq=""
	if ( /Score =/) (line_start=NR)

}

if ( /Score =/)
{	
	
	split($0, LINE_SPLIT, "=")
	split(LINE_SPLIT[3], LINE_SPLITTED, ",")
	while (sub(/ /, "", LINE_SPLITTED[1])){}
	evalue = LINE_SPLITTED[1]
	split(LINE_SPLIT[2], TEMP, " ")
	bits=TEMP[1]
	
	if (!/dentities/){getline}
	split($0, LINE_SPLIT, "%")
	
	while (sub(/[\(\)QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm=\:]/, "", LINE_SPLIT[1])){}
	split(LINE_SPLIT[1], LINE_SPLITTED, " ")
	identity = LINE_SPLITTED[length(LINE_SPLITTED)]


}



if (in_target_line){
    if (/^ +Length = [0-9]+ *$/) { target_length=$NF;     in_target_line=0    }
  else{
    if (FULL_TARGET){      subject=subject " " $0}
  }

}


if ( /^>/)   {
	if (FULL_TARGET){  subject = substr($0, 2)  }
	else{ 	           subject = substr($1, 2) }

	if (line_start==-1)
	{line_start = NR }
  in_target_line=1

}

  


}

END {
if (start_pos!=-1){

	if (ALI){
		COMMENT=COMMENT " " subj_seq " " query_seq 
		
	}
	if (SCORE)
	{	COMMENT=COMMENT " " evalue " " bits 
	}

	if (TLENGTH)
	{       COMMENT=COMMENT " " target_length
	}

	if (QLENGTH)
	{	COMMENT=COMMENT " " query_length 
	}


	if (REF_FILE)
	{	COMMENT=COMMENT " " ARGV[1] " " line_start " " (NR-1) 

	}

	strand="+"
	if (start_pos>end_pos){ temp=start_pos; start_pos=end_pos; end_pos=temp; strand="-"} #exchanging start with stop. using n_hits as a temp var

	if (! query){		query="UNKNOWN_QUERY" }

	if (FULL_QUERY){	gsub(/ +/, "%", query) }
	if (FULL_TARGET){	gsub(/ +/, "%", subject) }	
	
	print subject " " start_pos " " end_pos " " strand " "  query " " q_start_pos " " q_end_pos COMMENT


	}
}
