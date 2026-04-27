"report useful functions"

import os
import glob
import shutil
import re
#import configParser
from configparser import ConfigParser

def safecopy(file, path):
    "copy file to path"
    if not os.path.exists(file):
        raise Exception("ERROR! file not found: %s" % file)
    if not os.path.exists(path):
        os.makedirs(path)
    shutil.copy(file, path)


def loopcopy(pattern, path):
    '''
    pattern: shell regex pattern, "*.txt"
    path:    destination dir
    '''
    files = glob.glob(pattern)
    if files:
        for fig in glob.glob(pattern):
            safecopy(fig, path)
    else:
        raise Exception("ERROR! file not found: %s" % pattern)


def loopcopy5(pattern, path):#change_lzp
    '''
    pattern: shell regex pattern, "*.txt"
    path:    destination dir
    '''
    files = glob.glob(pattern)
    if files:
        for fig in glob.glob(pattern)[0:5]:#change_lzp
            safecopy(fig, path)
    else:
        raise Exception("ERROR! file not found: %s" % pattern)

def loopcopy20(pattern, path):#change_lzp
    '''
    pattern: shell regex pattern, "*.txt"
    path:    destination dir
    '''
    files = glob.glob(pattern)
    if files:
        for fig in glob.glob(pattern)[0:20]:#change_lzp
            safecopy(fig, path)
    else:
        raise Exception("ERROR! file not found: %s" % pattern)

def readtbl3(table, only5=False):
    "display xls table"
    content = []
    if not os.path.exists(table):
        raise Exception("ERROR! file not found: %s" % table)
    handle = open(table)
    handle.readline()
    if only5:
        for i in range(5):
            line = handle.readline()
            if line != "":
                content.append(line.strip().split('\t'))
    else:
        for line in handle:
            content.append(line.strip().split('\t'))
    handle.close()
    return content


def readtbl(table, only5=False):
    "display xls table"
    content = []
    if not os.path.exists(table):
        raise Exception("ERROR! file not found: %s" % table)
    handle = open(table)
    handle.readline()
    if only5:
        for i in range(5):
            line = handle.readline()
            if line != "":
                content2 = []
                for k in line.strip().split('\t'):
                    try:
                        kk = float(k)
                        pattern = '[0-9]*\\.[0-9]*'
                        if re.search(pattern, str(k)):
                            k=float(k)
                            float_i = format(k, '.3f')
                        else:
                            float_i = int(k)
                    except:
                        float_i = k
                    content2.append(float_i)
                content.append(content2)
                #content.append(line.strip().split('\t'))
    else:
        for line in handle:
            content2 = []
            for k in line.strip().split('\t'):
                try:
                    kk = float(k)
                    pattern = '[0-9]*\\.[0-9]*'
                    if re.search(pattern, str(k)):
                        k=float(k)
                        float_i = format(k, '.3f')
                    else:
                        float_i = int(k)
                except:
                    float_i = k
                content2.append(float_i)
            content.append(content2)
            #content.append(line.strip().split('\t'))
    handle.close()
    return content

def readtbl2(table, only5=False):
    if not os.path.exists(table):
        raise Exception("ERROR! file not found: %s" % table)
    handle = open(table)
    content = []
    header = handle.readline()
    header = header.replace('"','')
    header_line = ""
    for line in header.strip().split('\t'):
        header_line += "<th><b><yw>" + line + "</yw></b></th>\n"
    if only5:
        for i in range(5):
            line = handle.readline()
            line = line.replace('"','')
            line3 = ""
            if line != "":
                for line2 in line.strip().split('\t'):
                    try:
                        line5 = float(line2)
                        pattern = '[0-9]*\\.[0-9]*'
                        if re.search(pattern, str(line2)):
                            line2 = float(line2)
                            line4 = format(line2, '.3f')
                        else:
                            line4=int(line2)
                    except:
                        line4 = line2
                    line3 += "<td><yw>" + str(line4) + "</yw></td>\n"
                content.append(line3)
    else:
        for line in handle:
            line3 = ""
            for line2 in line.strip().split('\t'):
                try:
                    line5 = float(line2)
                    pattern = '[0-9]*\\.[0-9]*'
                    if re.search(pattern, str(line2)):
                        line2 = float(line2)
                        line4 = format(line2, '.3f')
                    else:
                        line4=int(line2)
                except:
                    line4 = line2
                line3 += "<td>" + str(line4) + "</td>\n"
            content.append(line3)
    handle.close()
    return header_line,content

def readtbl_go(table, only5=False):
    "display xls table"
    content = []
    if not os.path.exists(table):
        raise Exception("ERROR! file not found: %s" % table)
    handle = open(table)
    handle.readline()
    if only5:
        for i in range(5):
            line = handle.readline()
            if line != "":
                content2 = []
                for k in line.strip().split('\t'):
                    try:
                        kk = float(k)
                        pattern = '[0-9]*\\.[0-9]*'
                        if re.search(pattern, str(k)):
                            k=float(k)
                            float_i = format(k, '.3f')
                        else:
                            float_i = int(k)
                    except:
                        float_i = k
                    content2.append(float_i)
                line2 = content2[-1]
                line3 = re.sub(r"(\w*/\w*/\w*/\w*/\w*/\w*/)", r"\1 ",line2)
                content2[-1] = line3
                content.append(content2)
                #content.append(line.strip().split('\t'))
    else:
        for line in handle:
            content2 = []
            for k in line.strip().split('\t'):
                try:
                    kk = float(k)
                    pattern = '[0-9]*\\.[0-9]*'
                    if re.search(pattern, str(k)):
                        k=float(k)
                        float_i = format(k, '.3f')
                    else:
                        float_i = int(k)
                except:
                    float_i = k
                content2.append(float_i)
            line2 = content2[-1]
            line3 = re.sub(r"(\w*/\w*/\w*/\w*/\w*/\w*/)", r"\1 ",line2)
            content2[-1] = line3
            content.append(content2)
            #content.append(line.strip().split('\t'))
    handle.close()
    return content

def readraw(table, only5=False):
    "display xls table"
    content = []
    if not os.path.exists(table):
        raise Exception("ERROR! file not found: %s" % table)
    handle = open(table)
    if only5:
        for i in range(5):
            line = handle.readline()
            if line != "":
                content.append(line.strip().split('\t'))
    else:
        for line in handle:
            content.append(line.strip().split('\t'))
    handle.close()
    return content

def schfile(path, pattern="*.png"):
    "return file list"
    if not os.path.exists(path):
        raise Exception("ERROR! path not found: %s" % path)
    return [os.path.basename(fig) for fig in sorted(glob.glob(os.path.join(path, pattern)))]


def parse_cfg(cfg):
    "parse cfg file , return dict"
    if not os.path.exists(cfg):
        raise Exception("ERROR! config file not found: %s" % cfg)
    conf = ConfigParser.ConfigParser()
    conf.read(cfg)
    r = dict(conf.items('basic'))
    return r

