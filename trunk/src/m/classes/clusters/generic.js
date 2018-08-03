//GENERIC class definition
//
//   Usage:
//      generic=new generic();

function generic (){
	//properties 
	// {{{
	var args = Array.prototype.slice.call(arguments);
	var options = new pairoptions(args.slice(0,args.length));

	this.url=options.getfieldvalue('url','');
	this.np=options.getfieldvalue('np',3);
	this.codeversion=options.getfieldvalue('codeversion',20486);
	this.codepath=options.getfieldvalue('codepath','issmdir/bin');
	this.executionpath=options.getfieldvalue('executionpath','issmdir/execution');
	//}}}
	//methods
	this.disp= function(){// {{{
		console.log(sprintf('   generic class echo:'));
		console.log(sprintf('    url: "%s"',this.url));
		console.log(sprintf('    np: %i',this.np));
		console.log(sprintf('    codepath: "%s"',this.codepath));
		console.log(sprintf('    executionpath: "%s"',this.executionpath));
	}// }}}
	this.classname= function(){// {{{
		return "generic";
	}// }}}
	this.checkconsistency = function (md,solution,analyses) { //{{{
		if (cluster.np<1){
			md.checkmessage('number of processors should be at least 1');
		}
		if (isNaN(cluster.np)){
			md.checkmessage('number of processors should not be NaN!');
		}
	} //}}}
	this.BuildQueueScript = function (cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota) { // {{{

			//write queuing script 
			//what is the executable being called? 
			executable='issm.exe';

			fid=fopen(modelname+'.queue','w');
			fprintf(fid,'#!%s\n',cluster.shell);
			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s 2> %s.errlog >%s.outlog ',cluster.np,cluster.codepath,executable,EnumToString(solution),cluster.executionpath+'/'+dirname,modelname,modelname,modelname);					
			fclose(fid);
	} //}}}
	this.UploadAndRun = function (md,callbackfunction,callbackerrorfunction,callbackid,fid,toolkitsstring,solutionstring,name,runtimename) { //{{{
		if (!navigator.onLine) { //{{{
			$(callbackid).html(sprintf("%-16s", "NO CONNECTION")).prop("disabled", false);
			callbackerrorfunction();
			return;
		} //}}}
		var request = new XMLHttpRequest();
		$(callbackid).html(sprintf("%-16s", "CONNECTING...")).prop("disabled", true);
		request.position = 0; //Keep track of current parsing position in repsonseText
		request.timeout = 180000;
		request.ontimeout = function (event) { //{{{
			$(callbackid).html(sprintf("%-16s", "TIMEOUT")).prop("disabled", false);
			callbackerrorfunction();
		} //}}}
		request.onerror = function (event) { //{{{
			$(callbackid).html(sprintf("%-16s", "COULD NOT RUN")).prop("disabled", false);
			callbackerrorfunction();
		} //}}}
		request.upload.onprogress = function(event) { //{{{
			var progress = (event.loaded / event.total * 100).toFixed(0);
			$(callbackid).html(sprintf("%-20s", "UPLOADING: " + progress + "%"));
        } //}}}
		request.onprogress = function (event) { //{{{
			//Receive updates by parsing message length as a 32-bit hex string of form 0x*09ABCDEF))
			var startIndex = request.position;
			var endIndex = request.position + 10;
			if (request.responseText.length >= endIndex) { //Ensure entire hex string is loaded
				var chunkSize = parseInt(request.responseText.slice(startIndex, endIndex));
				startIndex = endIndex;
				endIndex = startIndex + chunkSize;
				if (chunkSize >= 1024) { //Arbitrary maximium size of message (Must be below minimium size of model results)
					$(callbackid).html(sprintf("%-20s", "DOWNLOADING: " + ((request.responseText.length - request.position) / chunkSize * 100).toFixed(0) + "%")).prop("disabled", true);
				}
				else if (request.responseText.length >= endIndex) { //Ensure entire chunk is loaded
					var responseChunk = request.responseText.slice(startIndex, endIndex);
					$(callbackid).html(responseChunk);
					request.position = endIndex;
				}
			}
		}; //}}}
		request.onload = function (event) { //{{{
			//get context to this.str2ab to avoid duplciation
			function str2ab(str) {
				var buf = new Uint8Array(str.length);
				for (var i=0, strLen=str.length; i < strLen; i++) {
					buf[i] = str.charCodeAt(i);
				}
				return buf;
			}
			var responseText = window.atob(request.responseText.slice(request.position + 10).replace(/\s/g, ''));
            var buffer = pako.inflate(str2ab(responseText));
			var returnBuffer = new Uint8Array(buffer);
			var returnBuffer_size = returnBuffer.byteLength;
			try {
				//Write result buffer to file for debugging. Filename and MIME type are optional.
				//writetofile(returnBuffer, "resultBuffer", "application/octet-stream");
				md.results = parseresultsfrombuffer(md,returnBuffer,returnBuffer_size);
				$(callbackid).html(sprintf("%-16s", "RUN")).prop("disabled", false);
				callbackfunction();
			}
			catch (e) {
				if (responseText.startsWith('Error')) {
					console.log(responseText);
					$(callbackid).html(sprintf("%-16s", "ISSM ERROR")).prop("disabled", false);
				}
				else {
					$(callbackid).html(sprintf("%-16s", "JS ERROR")).prop("disabled", false);
					console.log(e);
				}
				callbackerrorfunction();
			}
			
		}; //}}}
		
		var npbuffer = this.str2ab(md.cluster.np.toString());
		npbuffer = pako.deflate(npbuffer);
		var nplength = new Uint32Array(1);
		nplength[0] = npbuffer.byteLength;
		
		var codeversionbuffer = this.str2ab(md.cluster.codeversion.toString());
		codeversionbuffer = pako.deflate(codeversionbuffer);
		var codeversionlength = new Uint32Array(1);
		codeversionlength[0] = codeversionbuffer.byteLength;
		
		var runtimenamebuffer = this.str2ab(runtimename);
		runtimenamebuffer = pako.deflate(runtimenamebuffer);
		var runtimenamelength = new Uint32Array(1);
		runtimenamelength[0] = runtimenamebuffer.byteLength;
		
		var namebuffer = this.str2ab(name);
		namebuffer = pako.deflate(namebuffer);
		var namelength = new Uint32Array(1);
		namelength[0] = namebuffer.byteLength;
		
		var toolkitsbuffer = this.str2ab(toolkitsstring);
		toolkitsbuffer = pako.deflate(toolkitsbuffer);
		var toolkitslength = new Uint32Array(1);
		toolkitslength[0] = toolkitsbuffer.byteLength;
		
		var solutionbuffer = this.str2ab(solutionstring);
		solutionbuffer = pako.deflate(solutionbuffer);
		var solutionlength = new Uint32Array(1);
		solutionlength[0] = solutionbuffer.byteLength;
		
		var binbuffer = new Uint8Array(fid.rawbuffer()); //seems that 16 bits length could be incompatible.
		binbuffer = pako.deflate(binbuffer);
		var binlength = new Uint32Array(1);
		binlength[0] = binbuffer.byteLength;
		
		var data = new Blob([nplength,npbuffer,codeversionlength,codeversionbuffer,runtimenamelength,runtimenamebuffer,namelength,namebuffer,toolkitslength,toolkitsbuffer,solutionlength,solutionbuffer,binlength,binbuffer]);
		
		request.open("POST", this.url, true);
		request.responseType = 'application/octet-stream';
		request.send(data);
		
	} //}}}
	this.ab2str = function(buf) { //{{{
		return String.fromCharCode.apply(null, new Uint16Array(buf));
	} //}}}
	this.str2ab = function(str) { //{{{
		var buf = new Uint8Array(str.length);
		for (var i=0, strLen=str.length; i < strLen; i++) {
			buf[i] = str.charCodeAt(i);
		}
		return buf;
	} //}}}
}
