__all__ = ['entrez_eutils']

from time import sleep

# External Dependencies
import requests

def get_robust(url, count = 0, **kwargs):
    response = requests.get(url, **kwargs)
    try:
        response.raise_for_status()
    except Exception, e:        
        if count < 5:
            print("Error occured during HTTP Request (Error: %s), retry %d" % (str(e), count))
            sleep(10)
            return get_robust(url, count + 1, **kwargs)
        else:
            response.raise_for_status()
    return response

def post_robust(url, data, count = 0, **kwargs):
    response = requests.post(url, data = data, **kwargs)
    try:
        response.raise_for_status()
    except Exception, e:        
        if count < 5:
            print("Error occured during HTTP Request (Error: %s), retry %d" % (str(e), count))
            sleep(10)
            return get_robust(url, data, count + 1, **kwargs)
        else:
            response.raise_for_status()
    return response
