use crate::data::blastdb::ber::decode_integer;

const K_CLASS_MASK: u8 = 0xc0;
const K_CONSTRUCTED_MASK: u8 = 0x20;
const K_SHORT_TAG_MASK: u8 = 0x1f;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Class {
    Universal = 0,
    Application = 1,
    ContextSpecific = 2,
    Private = 3,
}

impl From<u8> for Class {
    fn from(value: u8) -> Self {
        match value {
            0 => Class::Universal,
            1 => Class::Application,
            2 => Class::ContextSpecific,
            3 => Class::Private,
            _ => unreachable!("ASN.1 class is stored in two bits"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum UniversalTag {
    Eoc = 0,
    Boolean = 1,
    Integer = 2,
    BitString = 3,
    OctetString = 4,
    Null = 5,
    ObjectIdentifier = 6,
    Utf8String = 12,
    Sequence = 16,
    Set = 17,
    PrintableString = 19,
    T61String = 20,
    Ia5String = 22,
    UtcTime = 23,
    GeneralizedTime = 24,
    BmpString = 30,
}

impl UniversalTag {
    pub fn from_tag_number(tag_number: u32) -> Option<Self> {
        match tag_number {
            0 => Some(Self::Eoc),
            1 => Some(Self::Boolean),
            2 => Some(Self::Integer),
            3 => Some(Self::BitString),
            4 => Some(Self::OctetString),
            5 => Some(Self::Null),
            6 => Some(Self::ObjectIdentifier),
            12 => Some(Self::Utf8String),
            16 => Some(Self::Sequence),
            17 => Some(Self::Set),
            19 => Some(Self::PrintableString),
            20 => Some(Self::T61String),
            22 => Some(Self::Ia5String),
            23 => Some(Self::UtcTime),
            24 => Some(Self::GeneralizedTime),
            30 => Some(Self::BmpString),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TagInfo {
    pub tag_class: Class,
    pub constructed: bool,
    pub tag_number: u32,
}

impl Default for TagInfo {
    fn default() -> Self {
        Self {
            tag_class: Class::Universal,
            constructed: false,
            tag_number: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Node {
    pub tag: TagInfo,
    pub value: Vec<u8>,
    pub children: Vec<Node>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodeError {
    msg: String,
}

impl DecodeError {
    pub fn new(msg: impl Into<String>) -> Self {
        Self { msg: msg.into() }
    }
}

impl std::fmt::Display for DecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.msg)
    }
}

impl std::error::Error for DecodeError {}

fn parse_tag(data: &[u8], length: usize, offset: &mut usize) -> Result<TagInfo, DecodeError> {
    if *offset >= length {
        return Err(DecodeError::new(
            "unexpected end of buffer while reading tag",
        ));
    }

    let first = data[*offset];
    *offset += 1;
    let mut info = TagInfo {
        tag_class: Class::from((first & K_CLASS_MASK) >> 6),
        constructed: (first & K_CONSTRUCTED_MASK) != 0,
        tag_number: 0,
    };

    let tag = first & K_SHORT_TAG_MASK;
    if tag != K_SHORT_TAG_MASK {
        info.tag_number = u32::from(tag);
        return Ok(info);
    }

    let mut continue_reading = true;
    let mut shift_count = 0;
    while continue_reading {
        if *offset >= length {
            return Err(DecodeError::new(
                "unexpected end of buffer while reading long tag",
            ));
        }
        let byte = data[*offset];
        *offset += 1;
        continue_reading = (byte & 0x80) != 0;
        info.tag_number = (info.tag_number << 7) | u32::from(byte & 0x7f);
        shift_count += 7;

        if shift_count > 28 {
            return Err(DecodeError::new("tag number is excessively large"));
        }
    }
    Ok(info)
}

fn parse_length(
    data: &[u8],
    length: usize,
    offset: &mut usize,
    indefinite: &mut bool,
) -> Result<usize, DecodeError> {
    if *offset >= length {
        return Err(DecodeError::new(
            "unexpected end of buffer while reading length",
        ));
    }
    let first = data[*offset];
    *offset += 1;
    if (first & 0x80) == 0 {
        *indefinite = false;
        return Ok(first as usize);
    }
    let count = first & 0x7f;
    if count == 0 {
        *indefinite = true;
        return Ok(0);
    }
    if usize::from(count) > std::mem::size_of::<usize>() {
        return Err(DecodeError::new("length uses more bytes than supported"));
    }
    if *offset + usize::from(count) > length {
        return Err(DecodeError::new(
            "unexpected end of buffer while reading long length",
        ));
    }
    let mut value = 0usize;
    for _ in 0..count {
        value = (value << 8) | usize::from(data[*offset]);
        *offset += 1;
    }
    *indefinite = false;
    Ok(value)
}

fn is_eoc(data: &[u8], length: usize, offset: usize) -> bool {
    offset + 1 < length && data[offset] == 0 && data[offset + 1] == 0
}

fn decode_impl(
    data: &[u8],
    length: usize,
    offset: &mut usize,
    stop_at_eoc: bool,
) -> Result<Vec<Node>, DecodeError> {
    let mut nodes = Vec::new();

    while *offset < length {
        if stop_at_eoc && is_eoc(data, length, *offset) {
            *offset += 2;
            break;
        }

        let tag = parse_tag(data, length, offset)?;

        let mut indefinite = false;
        let content_length = parse_length(data, length, offset, &mut indefinite)?;

        if !indefinite && *offset + content_length > length {
            return Err(DecodeError::new("content length exceeds available data"));
        }

        let mut node = Node {
            tag,
            value: Vec::new(),
            children: Vec::new(),
        };

        if tag.constructed {
            let end = if indefinite {
                length
            } else {
                *offset + content_length
            };
            node.children = decode_impl(data, end, offset, indefinite)?;

            if !indefinite && *offset != end {
                return Err(DecodeError::new(
                    "constructed element did not consume its content",
                ));
            }
        } else {
            let end = if indefinite {
                length
            } else {
                *offset + content_length
            };
            if indefinite {
                return Err(DecodeError::new(
                    "indefinite length used for primitive value",
                ));
            }
            node.value.extend_from_slice(&data[*offset..end]);
            *offset = end;
        }

        nodes.push(node);
    }

    Ok(nodes)
}

pub fn decode(data: &[u8]) -> Result<Vec<Node>, DecodeError> {
    let mut offset = 0;
    decode_impl(data, data.len(), &mut offset, false)
}

fn class_label(cls: Class) -> &'static str {
    match cls {
        Class::Universal => "Universal",
        Class::Application => "Application",
        Class::ContextSpecific => "Context-specific",
        Class::Private => "Private",
    }
}

fn universal_tag_name(tag_number: u32) -> &'static str {
    match UniversalTag::from_tag_number(tag_number) {
        Some(UniversalTag::Eoc) => "EOC",
        Some(UniversalTag::Boolean) => "BOOLEAN",
        Some(UniversalTag::Integer) => "INTEGER",
        Some(UniversalTag::BitString) => "BIT STRING",
        Some(UniversalTag::OctetString) => "OCTET STRING",
        Some(UniversalTag::Null) => "NULL",
        Some(UniversalTag::ObjectIdentifier) => "OBJECT IDENTIFIER",
        Some(UniversalTag::Utf8String) => "UTF8String",
        Some(UniversalTag::Sequence) => "SEQUENCE",
        Some(UniversalTag::Set) => "SET",
        Some(UniversalTag::PrintableString) => "PrintableString",
        Some(UniversalTag::T61String) => "T61String",
        Some(UniversalTag::Ia5String) => "IA5String",
        Some(UniversalTag::UtcTime) => "UTCTime",
        Some(UniversalTag::GeneralizedTime) => "GeneralizedTime",
        Some(UniversalTag::BmpString) => "BMPString",
        None => "",
    }
}

fn hex_dump(data: &[u8]) -> String {
    let mut out = String::new();
    for &b in data {
        if !out.is_empty() {
            out.push(' ');
        }
        out.push_str(&format!("{b:02x}"));
    }
    out
}

fn is_printable_ascii(data: &[u8]) -> bool {
    data.iter()
        .all(|&b| b == 0x0a || (0x20..=0x7e).contains(&b))
}

fn decode_oid(data: &[u8]) -> String {
    if data.is_empty() {
        return String::new();
    }

    let mut out = String::new();
    let first = data[0];
    out.push_str(&(first / 40).to_string());
    out.push('.');
    out.push_str(&(first % 40).to_string());

    let mut value = 0u32;
    let mut seen = false;
    for &byte in &data[1..] {
        value = (value << 7) | u32::from(byte & 0x7f);
        seen = true;
        if (byte & 0x80) == 0 {
            out.push('.');
            out.push_str(&value.to_string());
            value = 0;
            seen = false;
        }
    }

    if seen {
        return String::new();
    }

    out
}

fn decode_string(data: &[u8]) -> String {
    String::from_utf8_lossy(data).into_owned()
}

fn describe_value(node: &Node) -> String {
    if node.tag.constructed {
        return String::new();
    }

    if node.tag.tag_class != Class::Universal {
        return String::new();
    }

    match UniversalTag::from_tag_number(node.tag.tag_number) {
        Some(UniversalTag::Boolean) => {
            if node.value.len() == 1 {
                if node.value[0] != 0 {
                    "TRUE".to_string()
                } else {
                    "FALSE".to_string()
                }
            } else {
                String::new()
            }
        }
        Some(UniversalTag::Integer) => decode_integer(&node.value).to_string(),
        Some(UniversalTag::OctetString) => {
            if is_printable_ascii(&node.value) {
                format!("\"{}\"", decode_string(&node.value))
            } else {
                String::new()
            }
        }
        Some(UniversalTag::Null) => "NULL".to_string(),
        Some(UniversalTag::ObjectIdentifier) => decode_oid(&node.value),
        Some(UniversalTag::Utf8String)
        | Some(UniversalTag::PrintableString)
        | Some(UniversalTag::T61String)
        | Some(UniversalTag::Ia5String)
        | Some(UniversalTag::BmpString)
        | Some(UniversalTag::UtcTime)
        | Some(UniversalTag::GeneralizedTime) => decode_string(&node.value),
        Some(UniversalTag::BitString)
        | Some(UniversalTag::Sequence)
        | Some(UniversalTag::Set)
        | Some(UniversalTag::Eoc)
        | None => String::new(),
    }
}

pub fn print_node(node: &Node, os: &mut String, depth: usize) {
    let indent = "-".repeat(depth);

    os.push_str(&indent);
    os.push_str("Class: ");
    os.push_str(class_label(node.tag.tag_class));
    os.push_str(", Constructed: ");
    os.push_str(if node.tag.constructed {
        "true"
    } else {
        "false"
    });
    os.push_str(", Tag: ");
    os.push_str(&node.tag.tag_number.to_string());

    if node.tag.tag_class == Class::Universal {
        let name = universal_tag_name(node.tag.tag_number);
        if !name.is_empty() {
            os.push_str(" (");
            os.push_str(name);
            os.push(')');
        }
    }
    os.push('\n');

    if !node.value.is_empty() {
        let decoded = describe_value(node);
        if !decoded.is_empty() {
            os.push_str(&indent);
            os.push_str("  Decoded: ");
            os.push_str(&decoded);
            os.push('\n');
        } else {
            os.push_str(&indent);
            os.push_str("  String: ");
            os.push_str(&decode_string(&node.value));
            os.push('\n');
        }
        os.push_str(&indent);
        os.push_str("  Raw (");
        os.push_str(&node.value.len().to_string());
        os.push_str(" bytes): ");
        os.push_str(&hex_dump(&node.value));
        os.push('\n');
    }

    for child in &node.children {
        print_node(child, os, depth + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_sequence_with_integer_and_printable_string() {
        let nodes = decode(&[0x30, 0x08, 0x02, 0x01, 0x2a, 0x13, 0x03, b'a', b'b', b'c']).unwrap();
        assert_eq!(nodes.len(), 1);
        assert_eq!(nodes[0].tag.tag_class, Class::Universal);
        assert!(nodes[0].tag.constructed);
        assert_eq!(nodes[0].tag.tag_number, UniversalTag::Sequence as u32);
        assert_eq!(nodes[0].children.len(), 2);
        assert_eq!(nodes[0].children[0].value, vec![0x2a]);
        assert_eq!(nodes[0].children[1].value, b"abc");

        let mut out = String::new();
        print_node(&nodes[0], &mut out, 0);
        assert!(out.contains("Class: Universal, Constructed: true, Tag: 16 (SEQUENCE)\n"));
        assert!(out.contains("-Class: Universal, Constructed: false, Tag: 2 (INTEGER)\n"));
        assert!(out.contains("-  Decoded: 42\n"));
        assert!(out.contains("-  Decoded: abc\n"));
    }

    #[test]
    fn test_decode_long_tag_long_length_and_oid() {
        let nodes = decode(&[0x1f, 0x81, 0x00, 0x81, 0x03, 0x2a, 0x86, 0x48]).unwrap();
        assert_eq!(nodes[0].tag.tag_number, 128);
        assert_eq!(nodes[0].value, vec![0x2a, 0x86, 0x48]);
    }

    #[test]
    fn test_decode_indefinite_constructed_and_rejects_primitive() {
        let nodes = decode(&[0x30, 0x80, 0x01, 0x01, 0xff, 0x00, 0x00]).unwrap();
        assert_eq!(nodes[0].children.len(), 1);
        assert_eq!(
            nodes[0].children[0].tag.tag_number,
            UniversalTag::Boolean as u32
        );

        let err = decode(&[0x04, 0x80, b'a', 0x00, 0x00]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "indefinite length used for primitive value"
        );
    }

    #[test]
    fn test_decode_reports_cpp_error_messages() {
        assert_eq!(
            decode(&[0x30, 0x01, 0x02]).unwrap_err().to_string(),
            "unexpected end of buffer while reading length"
        );
        assert_eq!(
            decode(&[0x04, 0x02, 0x01]).unwrap_err().to_string(),
            "content length exceeds available data"
        );
        assert_eq!(
            decode(&[0x1f, 0x81, 0x82, 0x83, 0x84, 0x05])
                .unwrap_err()
                .to_string(),
            "tag number is excessively large"
        );
    }
}
