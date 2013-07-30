#!/bin/bash

valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes ./bin/dynamix 2> valgrind.tmp
