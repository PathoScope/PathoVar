def to_text(tag):
    return tag.get_text()

def to_text_strip(tag):
    return tag.get_text().strip()

def to_int(tag):
    return int(tag.get_text().strip())

def to_attr_value(tag):
    return tag.attrs['value']