#!/bin/bash

valgrind --leak-check=full --track-origins=yes --show-reachable=yes ./bin/dynamix 2> valgrind.tmp
