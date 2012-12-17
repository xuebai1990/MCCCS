def fileLen(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def readLastN(inputname,lastN):
    fin=open(inputname)
    #with open(inputname) as fin:
    step = 74
    while True:
        try:
            fin.seek(-(step*lastN), 2)
        except IOError:
            fin.seek(0)
        pos = fin.tell()
        lines = fin.read().splitlines()
        if len(lines) > lastN or pos == 0:
            if len(lines)<lastN:
                raise IOError
            else:
                return lines[-lastN]
        step*=2
    fin.close()

def tail( f, window=20 ):
    f.seek( 0, 2 )
    bytes= f.tell()
    size= window
    block= -1
    while size > 0 and bytes+block*1024  > 0:
        # If your OS is rude about small files, you need this check
        # If your OS does 'the right thing' then just f.seek( block*1024, 2 )
        # is sufficient
        if (bytes+block*1024 > 0):
            ##Seek back once more, if possible
            f.seek( block*1024, 2 )
        else:
            #Seek to the beginning
            f.seek(0, 0)
        data= f.read( 1024 )
        linesFound= data.count('\n')
        size -= linesFound
        block -= 1
    f.seek( block*1024, 2 )
    f.readline() # find a newline
    lastBlocks= list( f.readlines() )
    print lastBlocks[-window:]

class ReadFile:
    def __init__(self,inputname):
        self.__fin=file(inputname)
        self.__index=0
    def __del__(self):
        try:
            self.__fin.close()
        except AttributeError:
            pass
    def readNum(self):
        if self.__index == 0:
            line=self.__fin.readline()
            if len(line) == 0:
                raise EOFError
                return
            self.__fields=line.split();
            self.__index = len(self.__fields)
        self.__index = self.__index-1;
        self.__fields[-1-self.__index]=self.__fields[-1-self.__index].lower().replace('d','e')
        try:
            result=int(self.__fields[-1-self.__index])
        except ValueError:
            try:
                result=float(self.__fields[-1-self.__index])
            except ValueError:
                result="NAN"
        return result
    def readline(self,nline=1):
        if self.__index == 0:
            self.__fin.readline()
        else:
            self.__index = 0
