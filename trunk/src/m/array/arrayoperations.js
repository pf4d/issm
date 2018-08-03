function ArrayMax(array){ //{{{
	return Math.max.apply(null,array);
} //}}}
function ArrayMax2D(array){ //{{{
	
	var max=0;

	for (var i=0;i<array.length;i++){
		var subarray=array[i];
		max=Math.max(max,ArrayMax(subarray));
	}

	return max;
} //}}}
function ArrayMin(array){ //{{{
	return Math.min.apply(null,array);
} //}}}
function ArraySum(array){ //{{{
	var sum=0;
	for(var i=0;i<array.length;i++)sum+=array[i];
	return sum;
} //}}}
function ArrayXPY(){ //{{{
    if (arguments.length<2)throw Error("ArrayXPY error message: sum has to be for at least two arrays!");

	//check internal consistency of arrays provided!: 
	var firstarray=arguments[0];
	var firstsize=firstarray.length;
	
	for(var a=1;a<arguments.length;a++){
		var array=arguments[a];
		if(array.length!=firstsize)throw Error("ArrayXPY error message: arrays provided as arguments are not of the same length!");
	}

	//do the sum:
	var sum=NewArrayFill(firstsize,0);
	for(var a=0;a<arguments.length;a++){
		var array=arguments[a];
		for(var i=0;i<array.length;i++){
			sum[i]+=array[i];
		}
	}
	return sum;

} //}}}
function ArrayOr(){ //{{{
    if (arguments.length<2)throw Error("ArrayOr error message: sum has to be for at least two arrays!");

	//check internal consistency of arrays provided!: 
	var firstarray=arguments[0];
	var firstsize=firstarray.length;
	
	for(var a=1;a<arguments.length;a++){
		var array=arguments[a];
		if(array.length!=firstsize)throw Error("ArrayOr error message: arrays provided as arguments are not of the same length!");
	}

	//do the or:
	var or=NewArrayFill(firstsize,0);
	for(var a=0;a<arguments.length;a++){
		var array=arguments[a];
		for(var i=0;i<array.length;i++){
			or[i] = or[i] | array[i];
		}
	}
	return or;

} //}}}
function ArrayMin2D(array){ //{{{
	
	var min=ArrayMax2D(array);

	for (var i=0;i<array.length;i++){
		var subarray=array[i];
		min=Math.min(min,ArrayMin(subarray));
	}

	return min;
} //}}}
function ListToMatrix(list, elementsPerSubArray) { //{{{
	var matrix = [], i, k;

	for (i = 0, k = -1; i < list.length; i++) {
		if (i % elementsPerSubArray === 0) {
			k++;
			matrix[k] = [];
		}

		matrix[k].push(list[i]);
	}

	return matrix;
} //}}}
function MatrixToList(matrixin) { //{{{

	var matrix=matrixin;

	if (!IsArray(matrix[0])) return matrix;
	else{
		var width = matrix[0].length;
		var length = matrix.length;
		var list= new Array(width*length);

		for(var i=0;i<length;i++){
			for(var j=0;j<width;j++){
				list[i*width+j]=matrix[i][j];
			}
		}
		return list;
	}
} //}}}
function IsArray(object) { //{{{

	var type=Object.prototype.toString.call( object );
	if( type === '[object Array]' ) return 1;
	if( type === '[object Float64Array]' ) return 1;
	if( type === '[object Float32Array]' ) return 1;
	if( type === '[object Int32Array]' ) return 1;
	if( type === '[object Int16Array]' ) return 1;
	if( type === '[object Uint32Array]' ) return 1;
	if( type === '[object Uint16Array]' ) return 1;
	if( type === '[object Uint8Array]' ) return 1;
	return 0;

} //}}}
function ArrayNot(array) { //{{{

	var notarray=array;
	for (var i=0;i<array.length;i++)notarray[i]=-array[i];
	return notarray;
} //}}}
function ArrayCopy(array) { //{{{

	var copy=[];
	for(var i=0;i<array.length;i++)copy[i]=array[i];
	return copy;
} //}}}
function ArrayPow(array,coefficient) { //{{{

	var powarray=array;
	for (var i=0;i<array.length;i++)powarray[i]=Math.pow(array[i],coefficient);
	return powarray;
} //}}}
function ArraySqrt(array) { //{{{

	var sqrtarray=array;
	for (var i=0;i<array.length;i++)sqrtarray[i]=Math.sqrt(array[i]);
	return sqrtarray;
} //}}}
function ArrayScale(array,alpha) { //{{{

	for (var i=0;i<array.length;i++)array[i]=array[i]*alpha;

} //}}}
function ArrayMag(array1,array2) { //{{{

	var arraymag=NewArrayFill(array1.length,0);
	for (var i=0;i<array1.length;i++)arraymag[i]=Math.sqrt(Math.pow(array1[i],2)+Math.pow(array2[i],2));
	return arraymag;
} //}}}
function ArrayAnyNaN(array) { //{{{

    if(IsArray(array[0])){
        for(var i=0;i<array.length;i++){
            for(var j=0;j<array[0].length;j++){
                if (isNaN(array[i][j])) return 1;
            }
        }
    }
    else{
        for(var i=0;i<array.length;i++){
            if (isNaN(array[i])) return 1;
        }
    }
    return 0;
} //}}}
function ArrayUnique(arr) { //{{{

	return arr.reverse().filter(function (e, i, arr) {
		    return arr.indexOf(e, i+1) === -1;
	}).reverse();
} //}}}
function ArraySort(array) { //{{{

	return array.sort(function(a, b) {
		return a - b;
	});

} //}}}
function ArrayAnyEqual(array,value) { //{{{
	
	if(!isNaN(value)){
		for(var i=0;i<array.length;i++){
			if (array[i]==value)return 1;
		}
	}
	else{
		for(var i=0;i<array.length;i++){
			if (isNaN(array[i]))return 1;
		}
	}
	return 0;
} //}}}
function ArrayAnyBelowOrEqual(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]<=value)return 1;
	}
	return 0;
} //}}}
function ArrayAnyBelowStrict(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]<value)return 1;
	}
	return 0;
} //}}}
function ArrayAnyAboveOrEqual(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]>=value)return 1;
	}
	return 0;
} //}}}
function ArrayAnyAboveStrict(array,value) { //{{{

	for(var i=0;i<array.length;i++){
		if (array[i]>value)return 1;
	}
	return 0;
} //}}}
function ArrayAnd(array1,array2) { //{{{

	var array=array1;
	for (var i=0;i<array1.length;i++)array[i]=array1[i] & array2[i];
	return array;
} //}}}
function ArrayIsMember(array1,array2) { //{{{

	var array=NewArrayFill(array1.length,0);
	for (var i=0;i<array1.length;i++){
		for(var j=0;j<array2.length;j++){
			if (array1[i] == array2[j]){
				array[i]=1;
				break;
			}
		}
	}
	return array;
} //}}}
function NewArrayFill(size,value) { //{{{

	return new Array(size).fill(value);
} //}}}
function NewArrayFillIncrement(size,start,increment) { //{{{

	var array=new Array(size); 

	for(var i=0;i<size;i++){
		array[i]=start+i*increment;
	}

	return array;
} //}}}
function ArrayFind(array,value) { //{{{
	
	//find number of indices
	var count=0;
	for (var i=0;i<array.length;i++)if(array[i]==value)count++;

	//allocate:
	var indices= NewArrayFill(count,0);

	//fill in:
	count=0;
	for (var i=0;i<array.length;i++){
		if(array[i]==value){
			indices[count]=i;
			count++;
		}
	}
	return indices;
} //}}}
function ArrayFindNot(array,value) { //{{{
	
	//find number of indices
	var count=0;
	for (var i=0;i<array.length;i++)if(array[i]!=value)count++;

	//allocate:
	var indices= NewArrayFill(count,0);

	//fill in:
	count=0;
	for (var i=0;i<array.length;i++){
		if(array[i]!=value){
			indices[count]=i;
			count++;
		}
	}
	return indices;
} //}}}
function Create2DArray(rows,cols) { //{{{
	var arr = [];

	for (var i=0;i<rows;i++) {
		arr[i] = new Array(cols);
	}

	return arr;
} //}}}
function MapIsEmpty(map) { //{{{
	for (var key in map){
		if(map.hasOwnProperty(key)){
			return false;
		}
	}
	return true;
} //}}}
function clone(obj) {//{{{
	
	var copy;

	// Handle the 3 simple types, and null or undefined
	if (null == obj || "object" != typeof obj) return obj;

	// Handle Date
	if (obj instanceof Date) {
		copy = new Date();
		copy.setTime(obj.getTime());
		return copy;
	}

	// Handle Array
	if (obj instanceof Array) {
		copy = [];
		for (var i = 0, len = obj.length; i < len; i++) {
			copy[i] = clone(obj[i]);
		}
		return copy;
	}

	// Handle Object
	if (obj instanceof Object) {
		copy = {};
		for (var attr in obj) {
			if (obj.hasOwnProperty(attr)) copy[attr] = clone(obj[attr]);
		}
		return copy;
	}

	throw new Error("Unable to copy obj! Its type isn't supported.");
} //}}}
function FloatFix(pointer,size) {//{{{

	var buffer=new Float64Array(size);
	for(var i=0;i<size;i++)buffer[i]=pointer[i];
	return buffer;


} //}}}
function NullFix(pointer,value) {//{{{

	if(pointer==null)return value;
	else{
		//check that the pointer values are not null: 
		if(IsArray(pointer)){
			if(IsArray(pointer[0])){
				for(var i=0;i<pointer.length;i++){
					for(var j=0;j<pointer[0].length;j++){
						if(pointer[i][j]==null)pointer[i][j]=value;
					}
				}	
			}
			else{
				for(var i=0;i<pointer.length;i++){
					if(pointer[i]==null)pointer[i]=value;
				}
			}
		}
		return pointer;
	}

} //}}}
