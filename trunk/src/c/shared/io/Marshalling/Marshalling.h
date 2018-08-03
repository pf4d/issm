/*\file Marshalling.h
 *\brief: macros to help automate the marshalling, demarshalling, and marshalling size routines. 
 */

#ifndef _MARSHALLING_H_
#define _MARSHALLING_H_

enum marshall_directions{
	MARSHALLING_FORWARD,
	MARSHALLING_BACKWARD,
	MARSHALLING_SIZE
};

#define MARSHALLING_ENUM(EN)\
	int type_enum=EN;\
	if(marshall_direction==MARSHALLING_FORWARD){\
		memcpy(*pmarshalled_data,&type_enum,sizeof(int));\
		*pmarshalled_data+=sizeof(int);\
	}\
	else if(marshall_direction==MARSHALLING_SIZE){\
		*pmarshalled_data_size+=sizeof(int);\
	}\
	else if(marshall_direction==MARSHALLING_BACKWARD){\
		*pmarshalled_data+=sizeof(int);\
	}\
	else _error_("Wrong direction during the Marshall process");\


#define MARSHALLING(FIELD)\
	\
	if(marshall_direction==MARSHALLING_FORWARD){\
		memcpy(*pmarshalled_data,&FIELD,sizeof(FIELD));\
		*pmarshalled_data+=sizeof(FIELD);\
	}\
	else if(marshall_direction==MARSHALLING_SIZE){\
		*pmarshalled_data_size+=sizeof(FIELD);\
	}\
	else if(marshall_direction==MARSHALLING_BACKWARD){\
		memcpy(&FIELD,*pmarshalled_data,sizeof(FIELD));\
		*pmarshalled_data+=sizeof(FIELD);\
	}\
	else _error_("Wrong direction during the Marshall process");


#define MARSHALLING_ARRAY(FIELD,TYPE,SIZE) \
	\
	if(marshall_direction==MARSHALLING_FORWARD){\
		memcpy(*pmarshalled_data,FIELD,SIZE*sizeof(TYPE));\
		*pmarshalled_data+=SIZE*sizeof(TYPE);\
	}\
	else if(marshall_direction==MARSHALLING_SIZE){\
		*pmarshalled_data_size+=SIZE*sizeof(TYPE);\
	}\
	else if(marshall_direction==MARSHALLING_BACKWARD){\
		memcpy(FIELD,*pmarshalled_data,SIZE*sizeof(TYPE));\
		*pmarshalled_data+=SIZE*sizeof(TYPE);\
	}\
	else _error_("Wrong direction during the Marshall process");


#define MARSHALLING_DYNAMIC(FIELD,TYPE,SIZE) \
	\
	{\
		bool field_null=true;\
		if (FIELD)field_null=false;\
		MARSHALLING(field_null);\
		\
		if(!field_null){\
			if(marshall_direction==MARSHALLING_FORWARD){\
					memcpy(*pmarshalled_data,FIELD,SIZE*sizeof(TYPE));\
					*pmarshalled_data+=SIZE*sizeof(TYPE);\
			}\
			else if(marshall_direction==MARSHALLING_SIZE){\
				*pmarshalled_data_size+=SIZE*sizeof(TYPE);\
			}\
			else if(marshall_direction==MARSHALLING_BACKWARD){\
				FIELD=xNew<TYPE>(SIZE);\
				memcpy(FIELD,*pmarshalled_data,SIZE*sizeof(TYPE));\
				*pmarshalled_data+=SIZE*sizeof(TYPE);\
			}\
			else _error_("Wrong direction during the Marshall process");\
		}\
	}

#endif	
