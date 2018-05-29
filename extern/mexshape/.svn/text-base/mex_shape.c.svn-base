/*================================================================= *
 * MEX_SHAPE.C
 *     Gateway routine to interface with shapelib library.
 *
 *     This file evolved from shpdump.c as provided by the shapelib
 *     distribution.  Some flotsam and jetsam from that file are
 *     still present.
 *
 * The calling syntax is:
 *     [s, t] = mex_shapefile ( shapefile );
 *
 *     "s" is a structure array with at least three fields, "mx_data" and
 *     "my_data", which contain the vertices x and y data for each
 *     respective structure element, as well as the shape type.
 *
 *     "type" is a string defining the shape type.  This can be
 *     "Polygon"
 *
 * PARAMETERS:
 * Input:
 *   shapefile:
 *      character filename of "*.shp" file
 * Output:
 *   shape:
 *      structure array of shapefiles with their associated attributes.
 *   type:
 *      type of shapefile that was read.
 *
 * In case of an error, an exception is thrown.
 *
 *=================================================================*/
/* $Revision: 1.3 $ */
#include <string.h>
#include "shapefil.h"

#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
     
{ 

	/* 
	 * Construct error and warning messages using this buffer.  
	 * */
	char error_msg[500];



	int nShapeType, nEntities, i, j;


	int buflen;            /* length of input shapefile name */
	char *shapefile;       /* holder for input shapefile name */
	int status;            /* success or failure */

	/*
	 * Not used.
	 * */
	double adfMinBound[4], adfMaxBound[4];

	/*
	 * pointer to the shapefile
	 * */
	SHPHandle	hSHP;
	SHPObject	*psShape;



	/*
	 * handle for DBF file
	 * */
	DBFHandle     dbh;

	/*
	 * number of DBF attributes, records
	 * */
	int num_dbf_fields, num_dbf_records;


	/*
	 * This structure will hold a description of each field in the DBF file.
	 * */
	typedef struct DBF_Field_Descriptor {
		char          pszFieldName[12];
		DBFFieldType  field_type;
	} DBF_Field_Descriptor;


	DBF_Field_Descriptor *dbf_field;


	/*
	 * stores individual values from the DBF.
	 * */
	double dbf_double_val;
	int dbf_integer_val;
	char *dbf_char_val;

	char error_buffer[500];


	char **att_fnames;             /* holds name of fie */

	/* holds name of fields for the shape structure */
	char *shape_fnames[2] = { "x", "y" };      

	/* holds name of fields for the shape structure */
	char *outstruct_fnames[3]
		= { "Shape", "Attribute", "Type" };  


	/*
	 * Matlab structure to hold all output information.
	 * */
	mxArray *out_struct;

	/*
	 * Matlab structure that holds the point information.
	 *
	 * */
	mxArray *data_struct;
	
	/*
	 * Matlab structure that holds the attribute information.
	 * */
	mxArray *att_struct;
	
	/*
	 * 
	 * */
	mxArray *x_out, *y_out;

	/*
	 * temporary matlab array
	 * */
	mxArray *mxtmp;


	double *x_out_ptr, *y_out_ptr;

	/*
	 * Shortcuts into the output arrays.
	 * */
	double *mx_ptr, *my_ptr;

	int part_start;        /* start of a polygon part */
	int part_count;        /* how many vertices in a polygon part. */

	int dims[2];
	size_t sizebuf;

	/*
	 * This string will describe the type of shapefile.
	 * The possibilities are currently
	 *
	 *     MultiPoint
	 *     Point
	 *     Arc
	 *     Polygon
	 * */
	char *shapeTypeString;

	/*
	 * Initialize the dbf record.
	 * */
	dbf_field = NULL;



	/* Check for proper number of arguments */
	if (nrhs != 1) { 
		mexErrMsgTxt("One input argument is required."); 
	}
	if (nlhs != 1) { 
		mexErrMsgTxt("One output argument is required."); 
	}


	/*
	 * Make sure the input is a proper string.
	 * */
	if ( mxIsChar(prhs[0]) != 1 ) {
		mexErrMsgTxt("Shapefile parameter must be a string\n" );
	}
	if ( mxGetM(prhs[0]) != 1 ) {
		mexErrMsgTxt("Shapefile parameter must be a row vector, not a column string\n" );
	}

	buflen = mxGetN(prhs[0]) + 1; 
	shapefile = mxCalloc ( buflen, sizeof(char) );

	/*
	 * copy the string data from prhs[0] into a C string.
	 * */
	status = mxGetString ( prhs[0], shapefile, buflen );
	if ( status != 0 ) {
		mexErrMsgTxt ( "Not enough space for shapefile argument.\n" );
	}







	/* -------------------------------------------------------------------- */
	/*      Open the passed shapefile.                                      */
	/* -------------------------------------------------------------------- */
	hSHP = SHPOpen( shapefile, "rb" );
	if( hSHP == NULL ) {
		sprintf ( error_msg, "Unable to open:%s\n", shapefile );
		mexErrMsgTxt( error_msg );
	}




	/* -------------------------------------------------------------------- */
	/*      Get the needed information about the shapefile.                 */
	/* -------------------------------------------------------------------- */
	SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

	

	/*
	 * Make sure that we can handle the type.
	 * */
	switch ( nShapeType ) {
		case SHPT_MULTIPOINT:
		case SHPT_POINT:
		case SHPT_ARC:
		case SHPT_POLYGON:
			break;

		default:
			sprintf ( error_buffer, "Unhandled shape code %d (%s)\n", nShapeType, SHPTypeName ( nShapeType ) );
	        	mexErrMsgTxt( error_buffer ); 

	}
    
	/*
	 * Create the output shape type parameter.  
	 * */
	shapeTypeString = (char *) SHPTypeName ( nShapeType );


	/*
	 * Open the DBF in order to retrieve the number of fields and records.
	 * */
	dbh = DBFOpen (shapefile, "rb");
	num_dbf_fields = DBFGetFieldCount ( dbh );
	num_dbf_records = DBFGetRecordCount ( dbh );


	/*
	 * Allocate space for a description of each record, and populate it.
	 * I allocate space for two extra "dummy" records that go in positions
	 * 0 and 1.  These I reserve for the xy data.
	 * */
	dbf_field = (DBF_Field_Descriptor *) mxCalloc ( num_dbf_fields, sizeof ( DBF_Field_Descriptor ) );
	if ( dbf_field == NULL ) {
		mexErrMsgTxt("Memory allocation for DBF_Field_Descriptor failed."); 
	}


	for ( j = 0; j < num_dbf_fields; ++j ) {
		dbf_field[j].field_type = DBFGetFieldInfo ( dbh, j, dbf_field[j].pszFieldName, NULL, NULL );
	}




	/*
	 * Allocate space for the datapoint structure.
	 * */
	data_struct = mxCreateStructMatrix ( nEntities, 1, 2, (const char **)shape_fnames );


	/*
	 * Allocate space for the field names for the attributes.
	 * According to the API, each field name can have up to 12
	 * characters.
	 * */
	att_fnames = (char **) mxCalloc ( num_dbf_fields, sizeof(char *) );
	if ( att_fnames == NULL ) {
		mexErrMsgTxt("Memory allocation for attribute field names failed."); 
	}

	/*
	 * Copy the attribute names, create the matlab structure.
	 * */
	for ( j = 0; j < num_dbf_fields; ++j ) {
		att_fnames[j] = dbf_field[j].pszFieldName;
	}
	att_struct = mxCreateStructMatrix ( nEntities, 1, num_dbf_fields, (const char **)att_fnames );

    
	/* -------------------------------------------------------------------- */
	/*	Skim over the list of shapes, printing all the vertices.	*/
	/* -------------------------------------------------------------------- */
	for( i = 0; i < nEntities; i++ ) {

		psShape = SHPReadObject( hSHP, i );
		/*
		fprintf ( stdout, "Num parts = %d\n", psShape->nParts );
		*/


		/*
		 * Create the fields in this struct element.
		 * We will stick in one nan for each shape part.
		 */
		dims[0] = psShape->nVertices + psShape->nParts;
		dims[1] = 1;
		x_out = mxCreateNumericArray ( 2, dims, mxDOUBLE_CLASS, mxREAL );
		x_out_ptr = mxGetData ( x_out );
		y_out = mxCreateNumericArray ( 2, dims, mxDOUBLE_CLASS, mxREAL );
		y_out_ptr = mxGetData ( y_out );



		/*
		 * Just copy the verticies over.
		 * */
		sizebuf = mxGetElementSize ( x_out ) * psShape->nVertices;

		mx_ptr = x_out_ptr;
		my_ptr = y_out_ptr;


		for ( j = 0; j < psShape->nParts-1; ++j ) {

			part_count = psShape->panPartStart[j+1] - psShape->panPartStart[j];
			sizebuf = mxGetElementSize ( x_out ) * part_count;

			part_start = psShape->panPartStart[j];
			memcpy ( (void *) (mx_ptr), (void *) &(psShape->padfX[part_start]), sizebuf );
			memcpy ( (void *) (my_ptr), (void *) &(psShape->padfY[part_start]), sizebuf );

			/*
			 * Stick a nan here to separate the parts.
			 * */
			mx_ptr[part_count] = mxGetNaN();
			my_ptr[part_count] = mxGetNaN();

			/*
			 * Update the pointers to the next part.
			 * */
			mx_ptr += (part_count+1);
			my_ptr += (part_count+1);

		}


		/*
		 * Do the last one
		 *
		 * Special case if there is only a single point?  
		 * */
		if ( psShape->nParts == 0 ) {
			memcpy ( (void *) (mx_ptr), (void *) &(psShape->padfX[0]), sizebuf );
			memcpy ( (void *) (my_ptr), (void *) &(psShape->padfY[0]), sizebuf );
		} else {

			part_count = psShape->nVertices - psShape->panPartStart[psShape->nParts-1];
			sizebuf = mxGetElementSize ( x_out ) * part_count;
	
			part_start = psShape->panPartStart[psShape->nParts-1];
			memcpy ( (void *) (mx_ptr), (void *) &(psShape->padfX[part_start]), sizebuf );
			memcpy ( (void *) (my_ptr), (void *) &(psShape->padfY[part_start]), sizebuf );

			/*
			 * Stick a nan here to separate the parts.
			 * */
			mx_ptr[part_count] = mxGetNaN();
			my_ptr[part_count] = mxGetNaN();

		}


		mxSetField ( data_struct, i, "x", x_out );
		mxSetField ( data_struct, i, "y", y_out );






		/*
		 * Now do the attributes
		 * */
		for ( j = 0; j < num_dbf_fields; ++j ) {
			switch ( dbf_field[j].field_type ) {
				case FTString:
					dbf_char_val = (char *) DBFReadStringAttribute ( dbh, i, j );
					mxtmp = mxCreateString ( dbf_char_val );
					mxSetField ( att_struct, i, dbf_field[j].pszFieldName, mxtmp );
					break;

				case FTDouble:
					dbf_double_val = DBFReadDoubleAttribute ( dbh, i, j );
					mxtmp = mxCreateDoubleScalar ( dbf_double_val );
					mxSetField ( att_struct, i, dbf_field[j].pszFieldName, mxtmp );
					break;

				case FTInteger:
				case FTLogical:
					dbf_integer_val = DBFReadIntegerAttribute ( dbh, i, j );
					dbf_double_val = dbf_integer_val;
					mxtmp = mxCreateDoubleScalar ( dbf_double_val );
					mxSetField ( att_struct, i, dbf_field[j].pszFieldName, mxtmp );
					break;

				default:
					sprintf ( error_buffer, "Unhandled code %d, shape %d, record %d\n", dbf_field[j].field_type, i, j );
	        			mexErrMsgTxt("Unhandled code"); 

			}
		}


        SHPDestroyObject( psShape );

	}



	/*
	 * Clean up, close up shop.
	 * */
	SHPClose( hSHP );

	DBFClose ( dbh );




	/*
	 * Allocate space for the output structure.
	 * */
	out_struct = mxCreateStructMatrix ( 1, 1, 3, (const char **)outstruct_fnames );

	/*
	 * Set the fields properly.
	 * */
	mxSetField ( out_struct, 0, "Shape", data_struct );
	mxSetField ( out_struct, 0, "Attribute", att_struct );
	mxSetField ( out_struct, 0, "Type", mxCreateString ( shapeTypeString ) );

	plhs[0] = out_struct;
	return;

}










