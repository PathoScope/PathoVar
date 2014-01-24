__all__ = ['entrez_eutils']

# External Dependencies
import requests

def get_robust(url, count = 0, **kwargs):
    response = requests.get(url, **kwargs)
    try:
        response.raise_for_status()
    except Exception, e:        
        if count < 5:
            print("Error occured during HTTP Request (Error: %s), retry %d" % (e.text, count))
            sleep(10)
            return get_robust(url, count + 1, **kwargs)
        else:
            response.raise_for_status()
    return response

