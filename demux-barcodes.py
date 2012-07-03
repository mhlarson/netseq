#!/usr/bin/env python
import os,sys,string,argparse,shutil,multiprocessing,datetime,textwrap
"""
Demux barcodes using a configuration file with multiprocessing support.
3-prime barcode demuxing not yet implemented
"""
__author__ = "Calvin Jan"
__version__ = "0.2"
__date__ = "2011-09-05"

PROGNAME=__file__

def printer(inp):
    now = datetime.datetime.now()
    print >>sys.stdout, "%s [%s]: %s" % (PROGNAME, now, inp)
    return

def loadRead(inFile):
    '''returns a quartet of lines from the input file, or EOF if end of file is reached'''
    lines = []
    for i in range(4):
        line=inFile.readline()
        lines.append(line)
    if lines[3]=='': return "EOF"
    return lines

def demuxWorker(barcodeConfiguration,procnum,start,end,pipe):
    #start and end are either ints or tuples (tuples if index read)
    barcodeToCounts = dict()
    if barcodeConfiguration['barcodePosition']=='5prime':
        barcodeToOutputFile = dict()
        for barcode,filename in barcodeConfiguration['barcodeToReadFile'].iteritems():
            path_components=filename.split("/")[:-1]
            if not os.path.exists("/".join(path_components)): os.mkdir("/".join(path_components))
            barcodeToOutputFile[barcode]=open("%s._%s"%(filename,procnum),'w')
            barcodeToCounts[barcode]=0
        curfile=open(barcodeConfiguration['inputReadFile'])
        curfile.seek(start)
        while True:
            if curfile.tell()>=end: break
            lines = loadRead(curfile)
            if lines=='EOF': break
            try:
                barcodeToOutputFile[lines[1][:barcodeConfiguration['barcodeLength']]].write(string.join(lines,''))
                barcodeToCounts[lines[1][:barcodeConfiguration['barcodeLength']]]+=1
            except KeyError:
                barcodeToOutputFile['invalid'].write(string.join(lines,''))
                barcodeToCounts['invalid']+=1
        curfile.close()
        for barcode,filename in barcodeConfiguration['barcodeToReadFile'].iteritems():
            barcodeToOutputfile[barcode].close()
            x=open("%s._%s.complete"%(filename,procnum),'w')
            x.close()

    elif barcodeConfiguration['barcodePosition']=='index_read':
        barcodeToReadFile = dict()
        for barcode,filename in barcodeConfiguration['barcodeToReadFile'].iteritems():
            path_components=filename.split("/")[:-1]
            if not os.path.exists("/".join(path_components)): os.mkdir("/".join(path_components))
            barcodeToReadFile[barcode]=open("%s._%s"%(filename,procnum),'w')
            barcodeToCounts[barcode]=0
        barcodeToIndexFile = dict()
        for barcode,filename in barcodeConfiguration['barcodeToIndexFile'].iteritems():
            path_components=filename.split("/")[:-1]
            if not os.path.exists("/".join(path_components)): os.mkdir("/".join(path_components))
            barcodeToIndexFile[barcode]=open("%s._%s"%(filename,procnum),'w')
        read1file=open(barcodeConfiguration['inputReadFile'])
        read1file.seek(start[0])
        indexfile=open(barcodeConfiguration['inputIndexFile'])
        indexfile.seek(start[1])
        while True:
            if indexfile.tell()>=end[1]:break
            index_lines = loadRead(indexfile)
            read1_lines = loadRead(read1file)
            if index_lines=='EOF': break
            assert index_lines[0][:-2]==read1_lines[0][:-2]
            try:
                barcodeToReadFile[index_lines[1][:barcodeConfiguration['barcodeLength']]].write(string.join(read1_lines,''))
                barcodeToIndexFile[index_lines[1][:barcodeConfiguration['barcodeLength']]].write(string.join(index_lines,''))
                barcodeToCounts[index_lines[1][:barcodeConfiguration['barcodeLength']]]+=1
            except KeyError:
                barcodeToReadFile['invalid'].write(string.join(read1_lines,''))
                barcodeToIndexFile['invalid'].write(string.join(index_lines,''))
                barcodeToCounts['invalid']+=1
        indexfile.close()
        read1file.close()
        for barcode,filename in barcodeConfiguration['barcodeToReadFile'].iteritems():
            barcodeToReadFile[barcode].close()
            x=open("%s._%s.complete"%(filename,procnum),'w')
            x.close()
        for barcode,filename in barcodeConfiguration['barcodeToIndexFile'].iteritems():
            barcodeToIndexFile[barcode].close()
            x=open("%s._%s.complete"%(filename,procnum),'w')
            x.close()            
    pipe.send(barcodeToCounts)
    pipe.close()
    return

def estimate_lines(filename):
    f=open(filename)
    for i in range(4):
        f.readline()
    bytes_per_quartet = f.tell()
    f.seek(0,2)
    total_size = f.tell()
    f.close()
    approximate_lines = total_size/bytes_per_quartet
    approximate_lines *= 4
    return approximate_lines

def makeChunks(barcodeConfiguration, num_chunks):
    if barcodeConfiguration['barcodePosition']=='5prime' or barcodeConfiguration['barcodePosition']=='3prime':
        f = open(barcodeConfiguration['inputReadFile'])
        f.seek(0,2)
        file_size = f.tell()
        chunk_size = file_size/num_chunks
        real_chunks = [0]
        for i in range(1,num_chunks):
            f.seek(i*chunk_size)
            f.readline() # move to end of current line
            while True:
                last_position = f.tell()
                curline = f.readline()
                if curline[0]=='@': break
                else: pass
        real_chunks.append(last_position)
        assert len(real_chunks) == num_chunks
        ends = real_chunks[1:]
        ends.append(file_size)
        return zip(real_chunks,ends)
    else:
        approximate_lines=estimate_lines(barcodeConfiguration['inputReadFile'])
        lines_per_chunk=approximate_lines/num_chunks
        remainder = lines_per_chunk%4
        lines_per_chunk -= remainder
        read1file=open(barcodeConfiguration['inputReadFile'])
        indexfile=open(barcodeConfiguration['inputIndexFile'])
        read1file.seek(0,2)
        read1_size = read1file.tell()
        indexfile.seek(0,2)
        index_size = indexfile.tell()
        read1_chunk_size = read1_size/num_chunks
        read1_chunks = [0]
        index_chunks = [0]
        read1file.seek(0)
        indexfile.seek(0)
        for i in range(1,num_chunks):
            for n in xrange(lines_per_chunk):
                read1file.readline()
                indexfile.readline()
            read1_chunks.append(read1file.tell())
            index_chunks.append(indexfile.tell())
        read1_ends = read1_chunks[1:]
        read1_ends.append(read1_size)
        index_ends = index_chunks[1:]
        index_ends.append(index_size)
        starts = zip(read1_chunks,index_chunks)
        ends = zip(read1_ends,index_ends)
        return zip(starts,ends)


def main():
    printer("Running script: "+string.join(sys.argv, ' '))
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,epilog=textwrap.dedent("""\
                    Example barcodeConfig.txt
                    ----------------------------------------                                      
                        <Skip index file fields if barcode is integral to read>                 
                        # barcode location: 'index_read' # or '5prime' or '3prime'
                        barcodePosition = 'index_read'
                        # input sequence file                  
                        inputReadFile = '/home/user/s_1_1_sequence.txt'
                        # input index read file   
                        inputIndexFile = '/home/user/s_1_2_sequence.txt'
                        # barcode length
                        barcodeLength = 6 # length of barcode                                      
                        # list of barcodes
                        barcodes = ['ATCACG', 'TGACCA', 'GCCAAT', 'CAGATC']
                        # demuxing dictionary: specify where reads with each bc should be written
                        barcodeToReadFile =                                                      
                        {'invalid':'/home/user/invalid.txt',                                       
                        'ATCACG':'/home/user/ATCACG.txt',                                          
                        'TGACCA':'/home/user/TGACCA.txt',                                          
                        'GCCAAT':'/home/user/GCCAAT.txt',                                          
                        'CAGATC':'/home/user/CAGATC.txt'}
                        # demuxing dictionary: specify where index reads should be written
                        barcodeToIndexFile =                                                       
                        {'invalid':'/home/user/idx_invalid.txt',                                   
                        'ATCACG':'/home/user/idx_ATCACG.txt',                                      
                        'TGACCA':'/home/user/idx_TGACCA.txt',                                      
                        'GCCAAT':'/home/user/idx_GCCAAT.txt',                                      
                        'CAGATC':'/home/user/idx_CAGATC.txt'}"""))
    parser.add_argument("infile",type=str,metavar="barcodeConfig.txt",
                        help="Input file used to configure demux")
    parser.add_argument("-p","--processors",metavar="N",
                        help="Number of processors to use (default: 1)",
                        default=1)
    
    args    = parser.parse_args()

    barcodeFilename = args.infile
    num_cores = int(args.processors)
    printer("Loading demux options...")
    barcodeFile = open(barcodeFilename)
    barcodeContents = barcodeFile.read()
    barcodeConfiguration = dict()
    exec(barcodeContents) in barcodeConfiguration
    #NB: the barcode sequence will not be stripped.                                                                   
    assert barcodeConfiguration['barcodePosition']=='5prime' or barcodeConfiguration['barcodePosition']=='3prime' or barcodeConfiguration['barcodePosition']=='index_read'
    printer("Splitting jobs...")
    chunks = makeChunks(barcodeConfiguration,num_cores)
    processes = []
    pipes = []
    tempfiles=[]
    for barcode,filename in barcodeConfiguration['barcodeToReadFile'].iteritems():
        path_components=filename.split("/")[:-1]
        if not os.path.exists("/".join(path_components)): os.mkdir("/".join(path_components))
    if barcodeConfiguration['barcodePosition']=='index_read': 
        indextempfiles=[]
        for barcode,filename in barcodeConfiguration['barcodeToIndexFile'].iteritems():
            path_components=filename.split("/")[:-1]
            if not os.path.exists("/".join(path_components)): os.mkdir("/".join(path_components))
    barcodeToCounts={}
    counter = 1
    printer("Starting jobs...")
    for start,end in chunks:
        conA,conB=multiprocessing.Pipe()
        p=multiprocessing.Process(target=demuxWorker,args=(barcodeConfiguration,counter,start,end,conB))
        p.start()
        processes.append(p)
        pipes.append(conA)
        if  barcodeConfiguration['barcodePosition']=='5prime' or barcodeConfiguration['barcodePosition']=='3prime':
            tempfiles.extend(["%s._%s"%(fn,counter) for fn in barcodeConfiguration['barcodeToFile'].values()])
            for barcode in barcodeConfiguration['barcodeToFile'].keys(): barcodeToCounts[barcode]=0
        elif barcodeConfiguration['barcodePosition']=='index_read':
            tempfiles.extend(["%s._%s"%(fn,counter) for fn in barcodeConfiguration['barcodeToReadFile'].values()])
            indextempfiles.extend(["%s._%s"%(fn,counter) for fn in barcodeConfiguration['barcodeToIndexFile'].values()])
            for barcode in barcodeConfiguration['barcodeToReadFile'].keys(): barcodeToCounts[barcode]=0
        counter+=1
    for index,proc in enumerate(processes):
        proc.join()
        curcounts=pipes[index].recv()
        for barcode,count in curcounts.iteritems(): barcodeToCounts[barcode]=barcodeToCounts.get(barcode,0)+count
    for f in tempfiles:
        assert os.path.isfile("%s.complete"%f)

    printer('Done demuxing.  Catenating files...')
    if  barcodeConfiguration['barcodePosition']=='5prime' or barcodeConfiguration['barcodePosition']=='3prime':
        for barcode,filename in barcodeConfiguration['barcodeToFile'].iteritems():
            with open(filename, 'a') as destination:
                for n in range(num_procs):
                    shutil.copyfileobj(open("%s._%s"%(filename,n+1),'r'),destination)
    elif barcodeConfiguration['barcodePosition']=='index_read':
        for barcode,filename in barcodeConfiguration['barcodeToReadFile'].iteritems():
            with open(filename, 'a') as destination:
                for n in range(num_cores):
                    shutil.copyfileobj(open("%s._%s"%(filename,n+1),'r'),destination)
        for barcode,filename in barcodeConfiguration['barcodeToIndexFile'].iteritems():
            with open(filename, 'a') as destination:
                for n in range(num_cores):
                    shutil.copyfileobj(open("%s._%s"%(filename,n+1),'r'),destination)

    printer("Cleaning up temp files...")
    if  barcodeConfiguration['barcodePosition']=='5prime' or barcodeConfiguration['barcodePosition']=='3prime':
        for f in tempfiles:
            os.remove(f)
            os.remove("%s.complete"%f)
    elif barcodeConfiguration['barcodePosition']=='index_read':
        for f in tempfiles:
            os.remove(f)
            os.remove("%s.complete"%f)
        for f in indextempfiles:
            os.remove(f)
            os.remove("%s.complete"%f)

    printer("Barcode summary follows:")
    bcsummary=[]
    for bc,c in barcodeToCounts.iteritems():bcsummary.append((bc,c))
    bcsummary.sort(key=lambda x:x[1])
    bcsummary.reverse()
    for i in bcsummary: print "%s\t%s"%(i[0],i[1])
    printer("Done.")

if __name__=="__main__":
    main()
