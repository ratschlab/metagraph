#!/usr/bin/env python2.7
#
# coding: utf-8
#

import os
import sys
import re
import numpy as np
import random


__author__ = 'Mikhail Karasikov'

# awk '/^>/{print $0 ";longitude=16.3728" ";latitude=43.8608" ";city=Sarajevo"; next}{print}'


DIR = os.path.dirname(os.path.realpath(__file__))
cities_table = DIR + '/cities.tsv'


cities = np.loadtxt(cities_table, delimiter='\t', dtype=str)
columns = cities[0]
cities = cities[1:]


def dms2dd(direction, degrees, minutes=0, seconds=0):
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
    if direction == 'S' or direction == 'W':
        dd *= -1
    return dd;

def dd2dms(deg):
    d = int(deg)
    md = abs(deg - d) * 60
    m = int(md)
    sd = (md - m) * 60
    return [d, m, sd]

def parse_dms(dms):
    parts = re.split('[^\d\w]+', dms)

    lat_parts = parts[:len(parts) // 2]
    lng_parts = parts[len(parts) // 2:]

    lat = dms2dd(lat_parts[-1], *lat_parts[:-1])
    lng = dms2dd(lng_parts[-1], *lng_parts[:-1])
    return lat, lng

def generate_labels():
    labels = []

    try:
        location = random.choice(cities)

        for attribute, label in zip(columns, location):
            if attribute != 'Coordinates':
                labels.append('{}={}'.format(attribute, label))
            else:
                (lat, lng) = parse_dms(label)
                labels.append('latitude={}'.format(lat))
                labels.append('longitude={}'.format(lng))
    except:
        pass

    return labels


def main():
    if len(sys.argv) != 2:
        print("Usage: {} <input.fa>".format(sys.argv[0]))
        exit(1)

    with open(sys.argv[1]) as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                print(';'.join([line] + generate_labels()).replace(' ', ''))
            else:
                print(line)


if __name__ == '__main__':
    main()
