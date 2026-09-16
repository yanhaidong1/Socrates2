#!/usr/bin/env python

import argparse
import glob
import sys
import os


with open ('./test.txt','w') as opt:
    opt.write('test <- ' + f"{1}" + '\n' + 'test <- ' + '1')
