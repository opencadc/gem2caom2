#!/bin/bash

#instrument="GNIRS"
#for ii in GN-2011B-Q-7-193-010 GN-2011B-Q-63-101-006 GN-2011B-Q-68-116-013 GN-2011B-Q-63-126-005 GN-2011B-Q-7-193-044 GN-2011B-Q-63-126-005
#instrument="GRACES"
#for ii in GN-2015B-Q-1-12-1003 GN-CAL20150604-1000-1072 GN-CAL20150604-1000-1081 GN-CAL20150807-1000-1035
#instrument="NIFS"
#for ii in GN-2012B-Q-51-108-001 GN-2012B-Q-73-172-002 GN-2016B-Q-27-31-007 GN-2012A-Q-57-151-005 GN-2011B-Q-59-70-016 GN-2012A-Q-2-22-007
#instrument="GSAOI"
#for ii in GS-2013B-Q-61-8-008 GS-2012B-SV-499-21-002 GS-2013B-DD-1-13-002 GS-CAL20181023-5-001 GS-CAL20140109-3-009 GS-CAL20130201-3-017 GS-2013B-Q-26-19-004
#instrument="F2"
#for ii in GS-2014B-Q-60-8-071 GS-2015A-Q-60-126-001 GS-2014B-Q-17-69-002 GS-2014B-Q-17-53-030 GS-2014B-Q-60-11-012 GS-2017A-Q-28-280-012 GS-2016A-FT-18-54-002 GS-2016A-C-1-86-003 GS-2017B-Q-18-96-006 GS-2017B-Q-45-156-079 GS-2015A-Q-63-328-014 GS-2018B-SV-301-144-020 GS-2015A-Q-63-365-016 GS-2015A-Q-63-406-001
#instrument=GPI
#for ii in GS-2018A-FT-101-5-043 GS-CAL20150410-1-007 GS-CAL20140323-6-014 GS-CAL20160224-12-017 GS-2014A-SV-408-6-003 GS-CAL20140320-7-011 GS-CAL20140315-4-004 GS-CAL20140316-5-012
#instrument=Michelle
#for ii in GN-2008A-Q-43-6-002 GN-2009B-C-1-62-001 GN-CAL20060413-7-001 GN-CAL20080308-8-016 GN-2006A-Q-58-9-007 GN-2007A-C-11-234-001 GN-2010B-C-3-21-002 GN-2006A-C-14-49-002
#instrument=TReCS
#for ii in GS-CAL20050102-1-001 GS-2005A-Q-15-1-001 GS-2005B-Q-10-22-003 GS-2012A-Q-7-31-001 GS-2008A-C-5-35-002
#instrument=bHROS
#for ii in GS-2006B-Q-47-76-003 GS-2006B-Q-7-32-008 GS-2005B-SV-302-20-001 GS-2005B-SV-301-16-005
#instrument=hrwfs
#for ii in GS-2003A-Q-6-1-002 GS-CAL20030303-7-0006 GS-CAL20030730-10-006 GS-CAL20031218-1-034 GS-CAL20030105-2-0072
# instrument=Phoenix
# for ii in GS-2006B-C-8-2-0052 GS-2006A-DD-1-1-0073 GS-2003B-Q-51-27-0073 GS-2006A-C-10-1-0258
# instrument = OSCIR
#for ii in GS-2001B-Q-31-9-007
# instrument = Flamingos
# for ii in GS-2002B-DD-2-6-0230
# instrument = HOKUPAA
#for ii in GN-2002A-DD-1-319-591 GN-CAL20020424-1-007
# instrument = GMOS
for ii in GS-2003B-Q-5-29-003 GS-2005B-Q-27-4-001 GS-CAL20060125-1-002 GS-2005B-Q-27-33-001 GS-2008A-Q-41-26-013 GN2009AQ021-04 GS-2005B-Q-22-17-006 GS-2005B-Q-22-29-004 GS-2005B-Q-54-2-004 GS-2006A-Q-48-15-002 GS-2009A-Q-30-6-007 GN-2009B-Q-121-15-003 GN-2010A-Q-35-10-002 GN-2009B-Q-121-15-001 GN-2007B-Q-112-14-018
do
	echo $ii
        caom2-repo read --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo GEMINI $ii > "${ii}.in.xml"
        output=$(grep fits "${ii}.in.xml" | grep gemini | awk -F'/' '{print $2}' | awk -F'\.' '{print $1}')
	echo $output
        curl "https://archive.gemini.edu/jsonsummary/canonical/filepre=${output}" > "${ii}.json"
        curl "https://archive.gemini.edu/fullheader/${output}" > "${output}.fits.header"
done
