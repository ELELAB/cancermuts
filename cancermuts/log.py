# log.py - logging facilities for the cancermuts package
# (c) 2018 Matteo Tiberti <matteo.tiberti@gmail.com>
# This file is part of cancermuts
#
# cancermuts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Nome-Programma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>.


"""
metadata classes --- :mod:`cancermuts.log`
================================================================
Logging facilities for cancermuts

"""

import logging

logger_name = 'cancermuts'
date_format = "%d/%m/%Y,%H:%M"
message_format = "%(asctime)s %(levelname)s [%(name)s]> %(message)s"

def start_logging(logfile=None, level=logging.DEBUG, format=message_format, datefmt=date_format, silence_libs=True):
	if silence_libs:
		logging.getLogger("urllib3").setLevel(logging.WARNING)
		logging.getLogger("requests").setLevel(logging.WARNING)
	if logfile is not None:
		logging.basicConfig(filename=logfile, level=level, format=format, datefmt=datefmt)
	else:
		logging.basicConfig(                  level=level, format=format, datefmt=datefmt)

	logger = logging.getLogger(logger_name)
	return logger

def logger_init(function):
	def wrapper(*args, **kwargs):
		this_self = args[0]
		this_self.log = logging.getLogger('.'.join([logger_name, this_self.__class__.__name__]))
		function(*args, **kwargs)
	return wrapper
