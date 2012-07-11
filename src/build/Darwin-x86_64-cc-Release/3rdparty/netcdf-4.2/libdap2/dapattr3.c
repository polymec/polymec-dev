/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/libncdap3/dapattr3.c,v 1.14 2009/12/03 03:42:38 dmh Exp $
 *********************************************************************/

#include "ncdap3.h"

#define OCHECK(exp) if((ocstat = (exp))) goto done;

/* Forward */
static NCerror buildattribute(char*,nc_type,NClist*,NCattribute**);
static int mergedas1(NCDAPCOMMON*, OCconnection, CDFnode* dds, OCobject das);
static int isglobalname3(char* name);

#ifdef IGNORE
/* Extract attributes from the underlying oc objects
   and rematerialize them with the CDFnodes.
*/

NCerror
dapmerge3(NCDAPCOMMON* nccomm, CDFnode* node)
{
    unsigned int i;
    char* aname;
    unsigned int nvalues,nattrs;
    void* values;
    OCtype atype;
    OCerror ocstat = OC_NOERR;
    NCerror ncstat = NC_NOERR;
    NCattribute* att;

    if(node->dds == OCNULL) goto done;
    OCHECK(oc_inq_nattr(nccomm->conn,node->dds,&nattrs));
    if(nattrs == 0) goto done;
    if(node->attributes == NULL) node->attributes = nclistnew();
    for(i=0;i<nattrs;i++) {
	ocstat = oc_inq_attr(nccomm->conn,node->dds,i,
			           &aname,
			           &atype,
                                   &nvalues,
                                   &values);
	if(ocstat != OC_NOERR) continue; /* ignore */
        if(aname == NULL || nvalues == 0 || values == NULL)
	    continue; /* nothing to do */
	ncstat = buildattribute(aname,octypetonc(atype),
				      nvalues,values,&att);
	if(ncstat == NC_NOERR)
            nclistpush(node->attributes,(ncelem)att);
	nullfree(aname);
	oc_attr_reclaim(atype,nvalues,values);
    }
done:
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}

#endif

static NCerror
buildattribute(char* name, nc_type ptype,
               NClist* values, NCattribute** attp)
{
    NCerror ncstat = NC_NOERR;
    NCattribute* att;

    att = (NCattribute*)calloc(1,sizeof(NCattribute));
    MEMCHECK(att,NC_ENOMEM);
    att->name = nulldup(name);
    att->etype = ptype;

    att->values = values;

    if(attp) *attp = att;

    return THROW(ncstat);
}

#ifdef IGNORE
static NCerror
cvttype(nc_type etype, char** srcp, char** dstp)
{
    unsigned int typesize = nctypesizeof(etype);
    char* src = *srcp;
    char* dst = *dstp;

    switch (etype) {
    case NC_STRING: case NC_URL: {
	char* ssrc = *(char**)src;	
	*((char**)dst) = nulldup(ssrc);
	srcp += typesize;
	dstp += typesize;
    } break;
    
    default:
	if(typesize == 0) goto fail;
	memcpy((void*)dst,(void*)src,typesize);
	srcp += typesize;
	dstp += typesize;
	break;
    }
    return NC_NOERR;
fail:
    nclog(NCLOGERR,"cvttype bad value: %s",oc_typetostring(etype));
    return NC_EINVAL;
}
#endif

/*
Given a das attribute walk it to see if it
has at least 1 actual attribute (no recursion)
*/
static int
hasattribute3(OCconnection conn, OCobject dasnode)
{
    int i;
    OCerror ocstat = OC_NOERR;
    int tf = 0; /* assume false */
    unsigned int nsubnodes;
    OCtype ocsubtype;
    OCobject* subnodes = NULL;

    OCHECK(oc_inq_class(conn,dasnode,&ocsubtype));
    if(ocsubtype == OC_Attribute) return 1; /* this is an attribute */
    ASSERT((ocsubtype == OC_Attributeset));

    OCHECK(oc_inq_nsubnodes(conn,dasnode,&nsubnodes));
    OCHECK(oc_inq_subnodes(conn,dasnode,&subnodes));
    for(i=0;i<nsubnodes;i++) {
        OCobject subnode = subnodes[i];
        OCHECK(oc_inq_class(conn,subnode,&ocsubtype));
	if(ocsubtype == OC_Attribute) {tf=1; break;}
    }
done:
    nullfree(subnodes);
    return tf;
}

/*
Duplicate the oc merge das and dds code, but
modify to capture such things as "strlen" and "dimname".
*/

int
dapmerge3(NCDAPCOMMON* nccomm, CDFnode* ddsroot, OCobject dasroot)
{
    unsigned int i,j;
    NCerror ncerr = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    OCconnection conn = nccomm->oc.conn;
    unsigned int nsubnodes, nobjects;
    OCobject* dasobjects = NULL;
    NClist* dasglobals = nclistnew();
    NClist* dasnodes = nclistnew();
    NClist* dodsextra = nclistnew();
    NClist* varnodes = nclistnew();
    NClist* allddsnodes = ddsroot->tree->nodes;

    if(ddsroot == NULL || dasroot == NULL) return NC_NOERR;

    nobjects = oc_inq_nobjects(conn,dasroot);
    dasobjects = oc_inq_objects(conn,dasroot);
    
    /* 1. collect all the relevant DAS nodes;
          namely those that contain at least one
          attribute value.
          Simultaneously look for potential ambiguities
          if found; complain but continue: result are indeterminate.
          also collect globals and DODS_EXTRA separately.
    */
    for(i=0;i<nobjects;i++) {
	OCobject das = dasobjects[i];
	OCtype octype;
        char* ocname = NULL;
	int isglobal = 0;
	int hasattributes = 0;

        OCHECK(oc_inq_class(conn,das,&octype));
	if(octype == OC_Attribute) continue; /* ignore these for now*/

        OCHECK(oc_inq_name(conn,das,&ocname));
	OCHECK(oc_inq_nsubnodes(conn,das,&nsubnodes));

	isglobal = (ocname == NULL ? 0 : isglobalname3(ocname));

	/* catch DODS_EXTRA */
	if(isglobal && ocname != NULL && strcmp(ocname,"DODS_EXTRA")==0) {
	    nclistpush(dodsextra,(ncelem)das);
	    nullfree(ocname);
	    continue;
	}
	if(ocname == NULL || isglobal) {
            nclistpush(dasglobals,(ncelem)das);
	    nullfree(ocname);
	    continue;
	}
	hasattributes = hasattribute3(conn,das);
	if(hasattributes) {
	    /* Look for previously collected nodes with same name*/
            for(j=0;j<nclistlength(dasnodes);j++) {
	        OCobject das2 = (OCobject)nclistget(dasnodes,j);
		char* ocname2;
	        OCHECK(oc_inq_name(conn,das2,&ocname2));
		if(ocname2 == NULL || ocname == NULL) goto loop;
		if(strcmp(ocname2,"DODS")==0) goto loop;
	        if(strcmp(ocname,ocname2)==0)
		        nclog(NCLOGWARN,"nc_mergedas: potentially ambiguous DAS name: %s",ocname2);
loop:
		nullfree(ocname2);
	    }
	    nclistpush(dasnodes,(ncelem)das);
	}
	nullfree(ocname);
    }

    /* 2. collect all the leaf DDS nodes (of type NC_Primitive)*/
    for(i=0;i<nclistlength(allddsnodes);i++) {
	CDFnode* dds = (CDFnode*)nclistget(allddsnodes,i);
	if(dds->nctype == NC_Primitive) nclistpush(varnodes,(ncelem)dds);
    }

    /* 3. For each das node, lncate matching DDS node(s) and attach
          attributes to the DDS node(s).
          Match means:
          1. DAS->fullname :: DDS->fullname
          2. DAS->name :: DDS->fullname (support DAS names with embedded '.'
          3. DAS->name :: DDS->name
	  4. special case for DODS. Apply 1-3 on DODS parent.
    */
    for(i=0;i<nclistlength(dasnodes);i++) {
	OCobject das = (OCobject)nclistget(dasnodes,i);
	char* ocfullname = NULL;
	char* ocbasename = NULL;

	if(das == OCNULL) continue;
	OCHECK(oc_inq_name(conn,das,&ocbasename));
	if(strcmp(ocbasename,"DODS")==0) {
	    OCobject container;
   	    OCHECK(oc_inq_container(conn,das,&container));
	    if(container == OCNULL) {
	        ASSERT(container != OCNULL);
	    }
	    ocfullname = makeocpathstring3(conn,container,".");
	} else {
	    ocfullname = makeocpathstring3(conn,das,".");
	}
        for(j=0;j<nclistlength(varnodes);j++) {
	    CDFnode* dds = (CDFnode*)nclistget(varnodes,j);
	    char* ddsfullname = makecdfpathstring3(dds,".");
	    if(strcmp(ocfullname,ddsfullname)==0
	       || strcmp(ocbasename,ddsfullname)==0
	       || strcmp(ocbasename,dds->ocname)==0) {
		mergedas1(nccomm,conn,dds,das);
		/* remove from dasnodes list*/
		nclistset(dasnodes,i,(ncelem)NULL);
	    }
	    nullfree(ddsfullname);
	}
	nullfree(ocfullname);
	nullfree(ocbasename);
    }

    /* 4. Assign globals */
    for(i=0;i<nclistlength(dasglobals);i++) {
	OCobject das = (OCobject)nclistget(dasglobals,i);
	mergedas1(nccomm,conn,ddsroot,das);
    }

    /* 5. Assign DOD_EXTRA */
    for(i=0;i<nclistlength(dodsextra);i++) {
	OCobject das = (OCobject)nclistget(dodsextra,i);
	mergedas1(nccomm,conn,ddsroot,das);
    }

done: /* cleanup*/
    nullfree(dasobjects);
    nclistfree(dasglobals);
    nclistfree(dasnodes);
    nclistfree(dodsextra);
    nclistfree(varnodes);
    if(ocstat != OC_NOERR) ncerr = ocerrtoncerr(ocstat);
    return THROW(ncerr);
}

static int
mergedas1(NCDAPCOMMON* nccomm, OCconnection conn, CDFnode* dds, OCobject das)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    unsigned int i,j,k;
    unsigned int nsubnodes;
    OCobject* subnodes = NULL;
    OCobject* dodsnodes = NULL;
    unsigned int ndodsnodes;

    if(dds == NULL || das == OCNULL) return NC_NOERR; /* nothing to do */
    if(dds->attributes == NULL) dds->attributes = nclistnew();
    /* assign the simple attributes in the das set to this dds node*/
    OCHECK(oc_inq_nsubnodes(conn,das,&nsubnodes));
    OCHECK(oc_inq_subnodes(conn,das,&subnodes));
    for(i=0;i<nsubnodes;i++) {
	OCobject attnode = subnodes[i];
	OCtype octype, ocetype;
	char* ocname = NULL;
	unsigned int ocnvalues;
        OCHECK(oc_inq_name(conn,attnode,&ocname));	
        OCHECK(oc_inq_class(conn,attnode,&octype));
	if(octype == OC_Attribute) {
	    NCattribute* att = NULL;
	    NClist* stringvalues;
            OCHECK(oc_inq_primtype(conn,attnode,&ocetype));	
	    OCHECK(oc_inq_dasattr_nvalues(conn,attnode,&ocnvalues));
	    stringvalues = nclistnew();
	    for(j=0;j<ocnvalues;j++) {
		char* stringval;
	        OCHECK(oc_inq_dasattr(conn,attnode,j,&ocetype,&stringval));
	        nclistpush(stringvalues,(ncelem)stringval);
	    }
	    ncstat = buildattribute(ocname,
				    octypetonc(ocetype),
				    stringvalues,
				    &att);				
	    if(ncstat) goto done;
            nclistpush(dds->attributes,(ncelem)att);
	} else if(octype == OC_Attributeset
		  && (strcmp(ocname,"DODS")==0
		      || strcmp(ocname,"DODS_EXTRA")==0)) {
	    /* Turn the DODS special attributes into into
               special attributes for dds node */
	    OCHECK(oc_inq_nsubnodes(conn,attnode,&ndodsnodes));
	    OCHECK(oc_inq_subnodes(conn,attnode,&dodsnodes));
	    for(j=0;j<ndodsnodes;j++) {
		char* dodsname = NULL;
		char newname[4096];
	        OCobject attnode = dodsnodes[j];
	        NCattribute* att = NULL;
	        NClist* stringvalues;
	        OCHECK(oc_inq_class(conn,attnode,&octype));
		if(octype != OC_Attribute) continue;
                OCHECK(oc_inq_primtype(conn,attnode,&ocetype));	
	        OCHECK(oc_inq_dasattr_nvalues(conn,attnode,&ocnvalues));
	        stringvalues = nclistnew();
	        for(k=0;k<ocnvalues;k++) {
		    char* stringval;
	            OCHECK(oc_inq_dasattr(conn,attnode,k,&ocetype,&stringval));
	            nclistpush(stringvalues,(ncelem)stringval);
		}
	        OCHECK(oc_inq_name(conn,attnode,&dodsname));
		/* Compute new special name */
		strcpy(newname,"_DODS_");
		strcat(newname,dodsname);
	        ncstat = buildattribute(newname,
				        octypetonc(ocetype),
				        stringvalues,
				        &att);				
		if(ncstat) goto done;
		att->invisible = 1;
            	nclistpush(dds->attributes,(ncelem)att);

		/* Define extra semantics associated with DODS and DODS_EXTRA attribute */
		if(strcmp(dodsname,"strlen")==0) {
		    unsigned int maxstrlen = 0;
		    if(nclistlength(stringvalues) > 0) {
		        char* stringval = (char*)nclistget(stringvalues,0);
			if(0==sscanf(stringval,"%u",&maxstrlen)) maxstrlen = 0;
		    }
		    dds->dodsspecial.maxstrlen = maxstrlen;
#ifdef DEBUG
fprintf(stderr,"%s.maxstrlen=%d\n",dds->ocname,(int)dds->dodsspecial.maxstrlen);
#endif
		} else if(strcmp(dodsname,"dimName")==0) {
		    if(nclistlength(stringvalues) > 0) {
		        char* stringval = (char*)nclistget(stringvalues,0);
		        dds->dodsspecial.dimname = nulldup(stringval);
#ifdef DEBUG
fprintf(stderr,"%s.dimname=%s\n",dds->ocname,dds->dodsspecial.dimname);
#endif
		    } else dds->dodsspecial.dimname = NULL;
		} else if(strcmp(dodsname,"Unlimited_Dimension")==0) {
		    if(nccomm->cdf.recorddimname != NULL) {
		        nclog(NCLOGWARN,"Duplicate DODS_EXTRA:Unlimited_Dimension specifications");
		    } else if(nclistlength(stringvalues) > 0) {
		        char* stringval = (char*)nclistget(stringvalues,0);
			nccomm->cdf.recorddimname = nulldup(stringval);
#ifdef DEBUG
fprintf(stderr,"%s.Unlimited_Dimension=%s\n",dds->ocname,nccomm->cdf.recorddimname);
#endif
		    }
		} /* else ignore */
	        nullfree(dodsname);
	    }
	    nullfree(dodsnodes);
	}
        nullfree(ocname);
    }

done:
    nullfree(subnodes);
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}

#ifdef IGNORE
static NCerror
dodsextra3(NCDAPCOMMON* nccomm, NCattribute* att)
{
    int i,j;
    OCtype octype;
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;

    for(i=0;i<nclistlength(dodsextra);i++) {
 	OCobject das = (OCobject)nclistget(dodsextra,i);
	unsigned int ndodsnodes;
	OCobject* dodsnodes = NULL;
        OCHECK(oc_inq_class(conn,das,&octype));
	if(octype != OC_Attributeset) continue;
	/* Get the attributes within the DODS_EXTRA */
        OCHECK(oc_inq_nsubnodes(conn,das,&ndodsnodes));
        OCHECK(oc_inq_subnodes(conn,das,&dodsnodes));
        for(j=0;j<ndodsnodes;j++) {
	    OCobject extranode = dodsnodes[j];
   	    char* dodsname = NULL;
	    char* stringval;
	    unsigned int ocnvalues;
	    OCHECK(oc_inq_class(conn,extranode,&octype));
	    if(octype != OC_Attribute) continue;
	    OCHECK(oc_inq_name(conn,extranode,&dodsname));
	    OCHECK(oc_inq_dasattr_nvalues(conn,extranode,&ocnvalues));
	    if(strcmp(dodsname,"Unlimited_Dimension")==0 && ocnvalues > 0) {
	        OCHECK(oc_inq_dasattr(conn,extranode,0,NULL,&stringval));
		nccomm->cdf.recorddimname = stringval;
	    }
	    nullfree(dodsname);
	}
	nullfree(dodsnodes);
    }
done:
    return ncstat;
}
#endif

static int
isglobalname3(char* name)
{
    int len = strlen(name);
    int glen = strlen("global");
    char* p;
    if(len < glen) return 0;
    p = name + (len - glen);
    if(strcasecmp(p,"global") != 0)
	return 0;
    return 1;
}

