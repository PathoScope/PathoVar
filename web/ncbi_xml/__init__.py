ET = None
try:
    from lxml import etree as ET
except ImportError, e:
    try:
        from xml.etree import cElementTree as ET
    except ImportError, e:
        try:
            from xml.etree import ElementTree as ET
        except ImportError, e:
            print("Unable to import xml.etree.ElementTree, xml.etree.cElementTree, or lxml.etree. What gives?")
            raise e

def take_def(key, dict1, dict2):
    val1 = val2 = None
    fail = 0
    try:
        val1 = dict1[key]
    except KeyError, e:
        fail += 1
    try:
        val2 = dict2[key]
    except KeyError, e:
        fail += 2

    # Parity Check
    if fail == 1:
        return val2
    if fail == 2:
        return val2
    # Both failed. 
    if fail == 3:
        raise KeyError("take_def: %s Both values are not present!" % key)

    if(type(val1) == type(val2) == bool):
        return val1 or val2

    if val1 is not None:
        return val1
    elif val2 is not None:
        return val2
    
    return None

def try_tag(tag, func, default = None):
    try:
        return func(tag)
    except:
        return default

def to_text(tag):
    return ' '.join([child.text for child in tag.getiterator() if child is not None and child.text])

def to_text_strip(tag):
    try:
        return ' '.join([child.text.strip() for child in tag.getiterator() if child is not None and child.text and child.text.strip() != ''])
    except Exception, e:
        print(tag)
        print(tag.getchildren())
        raise e

def to_int(tag):
    return int(tag.text)

def to_attr_value(tag, key = None):
    if key is not None:
        return tag.get(key)
    else:
        return tag.items()
