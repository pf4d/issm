/*!\file Bucket.h
 * \brief: header file for Bucket object
 */

#ifndef _BUCKET_H
#define _BUCKET_H

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
#include "../../shared/io/Comm/IssmComm.h"
#include "../toolkitsenums.h"
/*}}}*/

/*how many ISSM_MPI_Isend requests does it take to transfer the contents of a bucket to another cpu?*/
#define MATRIXBUCKETSIZEOFREQUESTS 7 
#define VECTORBUCKETSIZEOFREQUESTS 5 
typedef enum {VECTOR_BUCKET, MATRIX_BUCKET} BucketType;
template <class doubletype> class Bucket: public Object{

	private: 
		int type; //either a VECTOR_BUCKET or MATRIX_BUCKET
		int m,n; /*size of local matrix we are storing*/
		/*row and column indices of the matrix we are storing*/
		int* idxm;
		int* idxn; 
		doubletype* values; /*local matrix*/
		InsMode mode; /*mode of insertion for this bucket*/

	public: 

		/*constructors, destructors: */
		Bucket(){ /*{{{*/
			this->Initialize();
		} /*}}}*/
		Bucket(int min,int* idxmin,int nin,int* idxnin,doubletype* valuesin,InsMode modein){ /*{{{*/

			this->Initialize();

			this->type=MATRIX_BUCKET;
			this->m=min;
			this->n=nin;
			this->mode=modein;
			if(this->m){
				this->idxm=xNew<int>(this->m); 
				xMemCpy(this->idxm,idxmin,this->m);
			}
			if(this->n){
				this->idxn=xNew<int>(this->n); 
				xMemCpy(this->idxn,idxnin,this->n);
			}
			if(this->m*this->n){
				this->values=xNew<doubletype>(this->n*this->m);
				xMemCpy(this->values,valuesin,this->n*this->m);
			}
		} /*}}}*/
		Bucket(int min,int* idxmin,doubletype* valuesin,InsMode modein){ /*{{{*/ 
			this->Initialize();

			this->type=VECTOR_BUCKET; 
			this->m=min;
			this->mode=modein;
			if(this->m){
				this->idxm=xNew<int>(this->m); 
				xMemCpy(this->idxm,idxmin,this->m);

				this->values=xNew<doubletype>(this->m);
				xMemCpy(this->values,valuesin,this->m);
			}
		} /*}}}*/
		~Bucket(){ /*{{{*/
			xDelete<int>(idxm);
			xDelete<int>(idxn);
			xDelete<doubletype>(values);
		} /*}}}*/
		void Initialize(void){ /*{{{*/

			this->type=0;
			this->m=0;
			this->n=0;
			this->idxm=NULL;
			this->idxn=NULL;
			this->values=NULL;
			mode=INS_VAL;

		} /*}}}*/

		/*object virtual functions definitions:*/
		void    Echo(){ /*{{{*/
			_printf_("Bucket echo (cpu #: "<<IssmComm::GetRank()<<")\n");
			_printf_("bucket type: " << type << "\n");
			_printf_("num rows: "<<this->m<<" num cols: "<<this->n << "\n");
		} /*}}}*/
		void    DeepEcho(){ /*{{{*/
			int i,j;

			_printf_("Bucket echo (cpu #: "<<IssmComm::GetRank()<<")\n");
			_printf_("bucket type: " << type << "\n");
			_printf_("num rows: "<<this->m<<" num cols: "<<this->n << "\n");
			if(type==MATRIX_BUCKET){
				for (i=0;i<this->m;i++){
					_printf_("row "<<this->idxm[i]<<", column indices: \n");
					for (j=0;j<this->n;j++){
						_printf_(" "<<this->idxn[j] << "\n");
					}
					_printf_("values: \n");
					for (j=0;j<this->n;j++){
						_printf_(" "<<this->values[m*i+j] << "\n");
					}
				}
			}
			else if(type==VECTOR_BUCKET){
				for (i=0;i<this->m;i++){
					_printf_("row "<<this->idxm[i]<<", value " << this->values[i] << "\n");
				}
			}
			else _error_("unknown type of bucket!");
		}
		/*}}}*/
		int     Id(){ /*{{{*/
			return -1;
		} /*}}}*/
		int     ObjectEnum(){ /*{{{*/
			return -1;
		} /*}}}*/
		Object *copy()        {/*{{{*/
			if (this->type==MATRIX_BUCKET) return new Bucket(this->m,this->idxm,this->n,this->idxn,this->values,this->mode);
			else if (this->type==VECTOR_BUCKET) return new Bucket(this->m,this->idxm,this->values,this->mode);
			else _error_("No Copy of Bucket because its type is invalid."); };
		/*}}}*/
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/

		/*specific routines of Bucket: */
		void SpawnBucketsPerCpu(DataSet* bucketsofcpu_i,int rank_i,int* rowranks){ /*{{{*/

			/*go through our idxm index of rows this bucket owns, and spawn buckets  
			 *if these rows belong to cpu rank_i. Use rowranks to determine this.*/
			for(int i=0;i<m;i++){
				if (rowranks[idxm[i]]==rank_i){
					/*This row belongs to cpu rank_i, so spawn a bucket with this row, and add it to the bucketsofcpu_i dataset: */
					if(type==MATRIX_BUCKET){
						bucketsofcpu_i->AddObject(new Bucket(1,idxm+i,n,idxn,values+n*i,mode));
					}
					else{
						bucketsofcpu_i->AddObject(new Bucket(1,idxm+i,values+i,mode));
					}
				}
			}

		}; /*}}}*/
		int BucketType(void){ /*{{{*/

			return type;
		};
		/*}}}*/
		void Marshall(int** prow_indices_forcpu,int** pcol_indices_forcpu,doubletype** pvalues_forcpu,int** pmodes_forcpu){ /*{{{*/

			/*intermediary: */
			int         i;
			int         j;

			/*buffers: */
			int        *row_indices_forcpu = NULL;
			int        *col_indices_forcpu = NULL;
			doubletype *values_forcpu      = NULL;
			int        *modes_forcpu       = NULL;

			/*initialize buffers: */
			row_indices_forcpu=*prow_indices_forcpu;
			col_indices_forcpu=*pcol_indices_forcpu;
			values_forcpu=*pvalues_forcpu;
			modes_forcpu=*pmodes_forcpu;

			/*fill buffers with out values and indices and modes: */
			for(i=0;i<m;i++){
				for(j=0;j<n;j++){
					row_indices_forcpu[i*n+j]=idxm[i];
					col_indices_forcpu[i*n+j]=idxn[j];
					values_forcpu[i*n+j]=values[i*n+j];
					modes_forcpu[i*n+j]=mode;
				}
			}

			/*increment buffer for next Bucket who will marshall his data: */
			row_indices_forcpu+=(m*n);
			col_indices_forcpu+=(m*n);
			values_forcpu+=(m*n);
			modes_forcpu+=(m*n);

			/*output modified buffers: */
			*prow_indices_forcpu=row_indices_forcpu;
			*pcol_indices_forcpu=col_indices_forcpu;
			*pvalues_forcpu=values_forcpu;
			*pmodes_forcpu=modes_forcpu;
		};
		/*}}}*/
		void Marshall(int** prow_indices_forcpu,doubletype** pvalues_forcpu,int** pmodes_forcpu){ /*{{{*/

			/*intermediary: */
			int         i;

			/*buffers: */
			int        *row_indices_forcpu = NULL;
			doubletype *values_forcpu      = NULL;
			int        *modes_forcpu       = NULL;

			/*initialize buffers: */
			row_indices_forcpu=*prow_indices_forcpu;
			values_forcpu=*pvalues_forcpu;
			modes_forcpu=*pmodes_forcpu;

			/*fill buffers with out values and indices and modes: */
			for(i=0;i<m;i++){
				row_indices_forcpu[i]=idxm[i];
				values_forcpu[i]=values[i];
				modes_forcpu[i]=mode;
			}

			/*increment buffer for next Bucket who will marshall his data: */
			row_indices_forcpu+=m;
			values_forcpu+=m;
			modes_forcpu+=m;

			/*output modified buffers: */
			*prow_indices_forcpu=row_indices_forcpu;
			*pvalues_forcpu=values_forcpu;
			*pmodes_forcpu=modes_forcpu;
		};
		/*}}}*/
		int MarshallSize(void){ /*{{{*/

			if(type==MATRIX_BUCKET){
				return m*n;
			}
			else{
				return m;
			}
		};
		/*}}}*/
};

#endif  /* _BUCKET_H */
