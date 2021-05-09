#include <cwp.h>
#include <su.h>
#include <segy.h>
#include <string.h>

char *sdoc[]={
" ",
" SUADDOBSERVATION - add observation system from sps/xps for 2D survey line",
" ",
" suaddobservation < stdin > stdout sht= rev= rel= [Optional parameters]",
" ",
" sht - sps shot points coordinates file",
" rev - sps receiver points coordinates file",
" rel - sps relational file",
" ",
" Optional paramters:",
" verbose=0  display sps information",
" ",
" Author:",
" Xiaobo.Foo from Chengdu University of Technology",
" Contact me: xiaobo_foo@126.com",
" ",
NULL};

segy tr; /* Define global variable for trace header */

typedef struct{
	char rId;		/* Record indentification */
	char name[16]; 	/* Line Name */
	char point[11];	/* Save ep and point index. Just for easy to read from file */
	unsigned int ep;/* Energy source point number, shot number */
	unsigned int pi;/* Point index */
	unsigned int pc;/* Point code */
	int stat;		/* Source static correction in milliseconds */
	float depth;	/* Source depth */
	int del; 		/* Datum elevation at source */
	int ut; 		/* uphole time*/
	float dep;   	/* Water depth at source */
	float x;		/* Source coordinate - X */
	float y;		/* Source coordinate - Y */
	float elev;		/* Surface elevation at source */
	int day;
	int hh;
	int mm;
	int ss;

}SPS;

typedef struct{
	char rId;		/* Record indentification */
	char tn[6];		/* Field tape number */
	char fnc[4];	/* string for field record number */
	int fn;			/* Field record number, int */
	char fic;		/* Field record increment, int */
	int fi;			/* Field record increment, int */
	char ic;		/* instrument code */
	char name[16];	/* Line name */
	char pnc[8];	/* Point number, same as SPS ep */
	int pn;			/* Point number, same as SPS ep */
	char pic;		/* Point index */
	int pi;			/* Point index */
	char c1c[4];	/* From channel */
	char c2c[4];	/* End channel */
	int c1;			/* From channel */
	int c2;			/* End channel */
	char cic;		/* Channel increment */
	int ci;			/* Channel increment */
	char name_[16];	/* Line name */
	char r1c[8];
	char r2c[8];
	int r1;			/* From receiver */
	int r2; 		/* End receiver */
	int ri;			/* Receiver index */

}XPS;

static int getNumShot(FILE *_FILE_);
static void getSourceInfo(FILE *_FILE_, SPS *sInfo, const int n);
static int getNumReceiver(FILE *_FILE_);
static void getReceiverInfo(FILE *_FILE_, SPS *sInfo, const int n);
static int getNumXPS(FILE *_FILE_);
static void getXPSInfo(FILE *_FILE_, XPS *sInfo, const int n);
static void printsps(SPS *sInfo, int n);
static void printxps(XPS *xInfo, int n);

int main(int argc, char *argv[])
{
	char *sht;  /* shot file name */
	char *rev;  /* receiver file name */
	char *rel;	/* relational file name */
	int nShot;  /* number of shot from sht file */
	int nRece;	/* number of receiver from the file */
	int nXPS;	/* number of XPS, should be equal to nShot */
	int verbose;
	SPS *sInfo = NULL;
	SPS *rInfo = NULL;
	XPS *xInfo = NULL;
	FILE *shtfp = NULL;
	FILE *revfp = NULL;
	FILE *relfp = NULL;	

	initargs(argc, argv);
	requestdoc(1);
	
	MUSTGETPARSTRING("sht",&sht);
	MUSTGETPARSTRING("rev",&rev);
	MUSTGETPARSTRING("rel",&rel);

	if(!getparint("verbose", &verbose)) verbose=0;

	checkpars();

	/* Read shot points coordinates */
	/* From the file, we can get the shot number */
	shtfp = efopen(sht,"r");
	nShot = getNumShot(shtfp);
	if(nShot>0)
	{
		sInfo = (SPS *)malloc(nShot*sizeof(SPS));
	}
	getSourceInfo(shtfp, sInfo, nShot);

	if(verbose) printsps(sInfo, nShot);
	efclose(shtfp);

	/* Read the relational file, and get the number of receiver for per shot */
	revfp = efopen(rev,"r");
	nRece = getNumReceiver(revfp);

	if(nRece>0)
	{
		rInfo = (SPS *)malloc(nRece*sizeof(SPS));
	}
	getReceiverInfo(revfp, rInfo, nRece);
	if(verbose)	printsps(rInfo, nRece);
	efclose(revfp);

	/* Read receiver points coordinates */
	relfp = efopen(rel, "r");
	nXPS = getNumXPS(relfp);
	if(nXPS>0)
	{
		xInfo = (XPS *)malloc(nXPS*sizeof(XPS));
	}

	getXPSInfo(relfp, xInfo, nXPS);

	if(verbose) printxps(xInfo, nXPS);

	efclose(relfp);
	
	/* Read seismic data and set coordinates */
	/* Get first trace */
	if (!gettr(&tr)) err("can't get first trace!");

	int minFldr = xInfo[0].fn;
	int minRece = xInfo[0].r1;

	do{
		/* do main process */
		/* according shot number to modify sx, sy, gx, gy, offset et al. */
		
		/* There is a bug, but I don't want fix it, now. Too tired. */	
		tr.scalel = 100;
		tr.scalco = 100; /* set scale for sx, sy, gx, gy coordinates*/
		tr.selev = tr.scalel * sInfo[tr.fldr-minFldr].elev; 
		tr.sx = tr.scalco * sInfo[tr.fldr-minFldr].x; 		/* source x coordinates */
		tr.sy = tr.scalco * sInfo[tr.fldr-minFldr].y; 		/* source y coordinates */

		int id = xInfo[tr.fldr-minFldr].r1-minRece+tr.tracf-1;
		tr.gelev = tr.scalel * rInfo[id].elev; 
		tr.gx = tr.scalco * rInfo[id].x; 		/* geophone x coordinates */
		tr.gy = tr.scalco * rInfo[id].y; 		/* geophone y coordinates */

		tr.offset = (int) sqrt(pow(tr.sx-tr.gx,2)+pow(tr.sy-tr.gy,2));

		puttr(&tr);
	}while(gettr(&tr));


	/* free memories */
	if(sInfo)
	{
		free(sInfo);
		sInfo = NULL;
	}

	if(rInfo)
	{
		free(rInfo);
		rInfo = NULL;
	}

	if(xInfo)
	{
		free(xInfo);
		xInfo = NULL;
	}

	return (CWP_Exit());
}


static int getNumShot(FILE *_FILE_)
{
	rewind(_FILE_); 
	char tmp;
	int nShot=0;
	int nLine = 0;
	do{
		fscanf(_FILE_,"%c",&tmp);
		if(tmp == '\n')
		{
			nLine++;
		}
	}while(!feof(_FILE_));

	nLine--;

	nShot = nLine - 20;

	rewind(_FILE_);

	return nShot;
}

/* There is a bug, but I can't change it. */
static void getSourceInfo(FILE *_FILE_, SPS *sInfo, const int n)
{
	rewind(_FILE_);

	char tmp;
	char tmpString[8];
	int nLine = 0;
	int i = 0;

	/* Skip file head */
	do{
		fscanf(_FILE_,"%c", &tmp);

		if(tmp=='\n')
		{
			nLine++;
		}

		if(nLine==20)
		{
			break;
		}

	}while(!feof(_FILE_));
	
	/* Read source information */

	for(i=0; i<n; i++)
	{
		/* Read record identification */
		fscanf(_FILE_, "%c", &sInfo[i].rId);
		if(sInfo[i].rId == '\n')
		{
			fscanf(_FILE_, "%c", &sInfo[i].rId);
		}

		/* Read Line name */
		fscanf(_FILE_, "%s", sInfo[i].name);

		/* Read point. Contained point number, index and code */
		fscanf(_FILE_, "%s", sInfo[i].point);

		/* Not convert to point number, index and code. Because, no use for me */

		/* Read Static correction */
		fscanf(_FILE_, "%d", &sInfo[i].stat);

		/* Read point depth */
		fscanf(_FILE_, "%f", &sInfo[i].depth);

		/* Read seismic datum */
		fscanf(_FILE_, "%d", &sInfo[i].del);

		/* Read uphole time */
		fscanf(_FILE_, "%d", &sInfo[i].ut);

		/* Read Water depth */
		fscanf(_FILE_, "%f", &sInfo[i].dep);

		/* Read coordinate */
		fscanf(_FILE_, "%f", &sInfo[i].x);
		fscanf(_FILE_, "%f", &sInfo[i].y);
		fscanf(_FILE_, "%f", &sInfo[i].elev);

		/* Read time */
		fscanf(_FILE_, "%d", &sInfo[i].day);
		fscanf(_FILE_, "%d", &sInfo[i].hh);
		fscanf(_FILE_, "%d", &sInfo[i].mm);
		fscanf(_FILE_, "%d", &sInfo[i].ss);
	}
	
	/* convert to point number, index and code. */
	for(i=0; i<n; i++)
	{
		nLine = strlen(sInfo[i].point);
		
		/* Convert to point number */
		strncpy(tmpString, sInfo[i].point, nLine-3);	/* Pick string */
		sInfo[i].ep = atoi(tmpString); /* String to int */
	
		/* Convert to point index */
		sInfo[i].pi = atoi(&sInfo[i].point[nLine-3]); /*atoi(tmpString);*/
		sInfo[i].pc = atoi(&sInfo[i].point[nLine-1]);/*atoi(tmpString);*/
	}

}


static int getNumReceiver(FILE *_FILE_)
{

	return getNumShot(_FILE_);
}

static void getReceiverInfo(FILE *_FILE_, SPS *sInfo, const int n)
{
	getSourceInfo(_FILE_, sInfo, n);
}

static int getNumXPS(FILE *_FILE_)
{
	return getNumShot(_FILE_);
}

static void getXPSInfo(FILE *_FILE_, XPS *sInfo, const int n)
{
	rewind(_FILE_);

	char tmp;
	int nLine = 0;
	int i = 0;
	int j = 0;
	int counter=0;

	/* Skip file head */
	do{
		fscanf(_FILE_,"%c", &tmp);

		if(tmp=='\n')
		{
			nLine++;
		}

		if(nLine==20)
		{
			break;
		}

	}while(!feof(_FILE_));
	
	/* Read XPS information */
	for(i=0; i<n; i++)
	{
		/* Read record identification */
		fscanf(_FILE_, "%c", &sInfo[i].rId);
		if(sInfo[i].rId == '\n')
		{
			fscanf(_FILE_, "%c", &sInfo[i].rId);
		}
		
		/* Read Field tape number */
		counter=0;
		for(j=0; j<6; j++)
		{
			fscanf(_FILE_,"%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].tn[counter++]=tmp;
			}
		}
		sInfo[i].tn[counter] = '\0';

		/* Read Field record information */
		counter=0;
		for(j=0; j<4; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].fnc[counter++] = tmp;
			}
		}
		sInfo[i].fnc[counter] = '\0';
		/* Convert to int */
		sInfo[i].fn = atoi(sInfo[i].fnc);

		fscanf(_FILE_, "%c", &sInfo[i].fic);

		sInfo[i].fi = atoi(&sInfo[i].fic); /* String to int */
		/* Instrument code */
		fscanf(_FILE_, "%c", &sInfo[i].ic);
		
		/* Read name 1*/
		counter=0;
		for(j=0; j<16; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].name[counter++] = tmp;
			}
		}
		sInfo[i].name[counter] = '\0';

		/* Read point number */
		counter=0;
		for(j=0; j<8; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].pnc[counter++] = tmp;
			}
		}
		sInfo[i].pnc[counter] = '\0';
		sInfo[i].pn = atoi(sInfo[i].pnc);
		

		fscanf(_FILE_,"%c", &sInfo[i].pic);
		sInfo[i].pi = atoi(&sInfo[i].pic);

		/* Read start channel */
		counter=0;
		for(j=0; j<4; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].c1c[counter++] = tmp;
			}
		}
		sInfo[i].c1c[counter] = '\0';

		sInfo[i].c1 = atoi(sInfo[i].c1c);

		/* Read end channel */
		counter=0;
		for(j=0; j<4; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].c2c[counter++] = tmp;
			}
		}
		sInfo[i].c2c[counter] = '\0';

		sInfo[i].c2 = atoi(sInfo[i].c2c);

		fscanf(_FILE_,"%c", &sInfo[i].cic);
		sInfo[i].ci = atoi(&sInfo[i].cic);

		/* Read name 1*/
		counter=0;
		for(j=0; j<16; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].name_[counter++] = tmp;
			}
		}
		sInfo[i].name_[counter] = '\0';

		/* Read start receiver */
		counter=0;
		for(j=0; j<8; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].r1c[counter++] = tmp;
			}
		}
		sInfo[i].r1c[counter] = '\0';
		sInfo[i].r1 = atoi(sInfo[i].r1c);

		/* Read start receiver */
		counter=0;
		for(j=0; j<8; j++)
		{
			fscanf(_FILE_, "%c", &tmp);
			if(tmp!=' ')
			{
				sInfo[i].r2c[counter++] = tmp;
			}
		}
		sInfo[i].r2c[counter] = '\0';
		sInfo[i].r2 = atoi(sInfo[i].r2c);

		fscanf(_FILE_, "%d", &sInfo[i].ri);
	}
}

static void printsps(SPS *sInfo, int n)
{
	int i;

	for (i=0; i<n; i++)
	{
		fprintf(stderr,"%d: %c%s\t%d%dV%d\t",i+1,sInfo[i].rId,sInfo[i].name,sInfo[i].ep,sInfo[i].pc,sInfo[i].pc);
		fprintf(stderr,"%d\t%.1f\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%d\t%d\t%d\t%d\n",\
						sInfo[i].stat, sInfo[i].depth, sInfo[i].del,\
						sInfo[i].ut,   sInfo[i].dep,   sInfo[i].x, \
						sInfo[i].y,    sInfo[i].elev,  sInfo[i].day,\
						sInfo[i].hh,    sInfo[i].mm,     sInfo[i].ss);
	}
}

static void printxps(XPS *xInfo, int n)
{
	int i;

	for(i=0; i<n; i++)
	{
		fprintf(stderr,"%d: %c%s\t",i+1, xInfo[i].rId, xInfo[i].tn);
		fprintf(stderr,"%d%d%c%s\t",xInfo[i].fn, xInfo[i].fi, xInfo[i].ic, xInfo[i].name);
		fprintf(stderr,"%d%d\t",xInfo[i].pn, xInfo[i].pi);
		fprintf(stderr,"%d\t%d%d%s\t",xInfo[i].c1, xInfo[i].c2, xInfo[i].ci, xInfo[i].name_);
		fprintf(stderr,"%d\t%d%d\n",xInfo[i].r1,xInfo[i].r2,xInfo[i].ri);
	}
}

















