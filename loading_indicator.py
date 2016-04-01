import sys

#MARK: Console Loading Indicator

class ConsoleColor:
	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

_total_loading_amount = 0
_loading_progress = 0
_has_loading_indicator = False
def _update_loading_indicator():
	global _has_loading_indicator, _total_loading_amount, _loading_progress
	text = ""
	if _total_loading_amount == 0:
		text = "\r"
		if _has_loading_indicator:
			sys.stdout.write(text)
			sys.stdout.flush()
	else:
		text = ConsoleColor.darkcyan + "{:.0f}% ".format(float(_loading_progress) / _total_loading_amount * 100.0) + ConsoleColor.end
		if not _has_loading_indicator:
			text = "\r" + text
		sys.stdout.write(text)
		sys.stdout.flush()


def add_loading_data(amt_load):
	global _total_loading_amount
	_total_loading_amount += amt_load
	_update_loading_indicator()

def update_progress(amt_loaded):
	global _loading_progress
	_loading_progress += amt_loaded
	_update_loading_indicator()

def clear_loading_data():
	global _total_loading_amount, _loading_progress
	_total_loading_amount = 0
	_loading_progress = 0
	_update_loading_indicator()