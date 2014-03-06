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


def to_text(tag):
    return ''.join([child.text for child in tag.getiterator()])

def to_text_strip(tag):
    try:
        return ''.join([child.text.strip() for child in tag.getiterator() if child.text and child.text.strip() != ''])
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
