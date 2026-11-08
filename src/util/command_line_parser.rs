use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::util::options::OptionValue;

fn raw_error(message: Vec<u8>) -> std::io::Error {
    crate::util::io::RawDiagnosticError(message).into_io()
}
fn source_os_bytes(value: &std::ffi::OsStr) -> std::io::Result<std::borrow::Cow<'_, [u8]>> {
    crate::util::io::source_os_bytes(value).map_err(std::io::Error::other)
}
fn source_bytes_os(value: &[u8]) -> std::io::Result<std::ffi::OsString> {
    crate::util::io::source_bytes_os(value).map_err(std::io::Error::other)
}
#[derive(Clone)]
enum Argument {
    Bytes(Vec<u8>),
    Os(std::ffi::OsString),
}
impl Argument {
    fn bytes(&self) -> std::io::Result<std::borrow::Cow<'_, [u8]>> {
        match self {
            Self::Bytes(s) => Ok(std::borrow::Cow::Borrowed(s)),
            Self::Os(s) => source_os_bytes(s),
        }
    }
    fn attached(&self) -> std::io::Result<Self> {
        match self {
            Self::Bytes(s) => Ok(Self::Bytes(s[2..].to_vec())),
            Self::Os(s) => {
                // Only split after the two ASCII option-prefix code units.
                // Both are complete boundaries in OsStr's encoded representation.
                #[cfg(unix)]
                {
                    use std::os::unix::ffi::{OsStrExt, OsStringExt};
                    Ok(Self::Os(std::ffi::OsString::from_vec(
                        s.as_bytes()[2..].to_vec(),
                    )))
                }
                #[cfg(windows)]
                {
                    use std::os::windows::ffi::{OsStrExt, OsStringExt};
                    Ok(Self::Os(std::ffi::OsString::from_wide(
                        &s.encode_wide().skip(2).collect::<Vec<_>>(),
                    )))
                }
                #[cfg(not(any(unix, windows)))]
                {
                    let encoded = s.as_encoded_bytes();
                    let tail =
                        unsafe { std::ffi::OsStr::from_encoded_bytes_unchecked(&encoded[2..]) };
                    Ok(Self::Os(tail.to_owned()))
                }
            }
        }
    }
}

/// Typed option readers. Raw entrypoints preserve Source argv bytes; OS readers
/// additionally allow path cells to retain their native representation.
pub trait ReadOption: Clone + 'static {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String>;
    fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
        let values = v
            .iter()
            .map(|s| String::from_utf8(s.clone()).map_err(std::io::Error::other))
            .collect::<std::io::Result<Vec<_>>>()?;
        Self::read_option(dst, &values).map_err(|message| raw_error(message.into_bytes()))
    }
    fn read_option_os(dst: &mut Self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
        let values = v
            .iter()
            .map(|s| source_os_bytes(s).map(|v| v.into_owned()))
            .collect::<std::io::Result<Vec<_>>>()?;
        Self::read_option_bytes(dst, &values)
    }
    fn check_pcount(v: &[String], min_count: i32) -> bool {
        let _ = min_count;
        v.len() == 1
    }
    fn check_pcount_len(count: usize, min_count: i32) -> bool {
        Self::check_pcount(&vec![String::new(); count], min_count)
    }
    fn check_present(&self) -> bool {
        false
    }
    fn set_option_default(dst: &mut Self, value: &Self) {
        *dst = value.clone();
    }
    fn set_base_ptr(&mut self, _base: &OptionBase) {}
}

fn c_input(value: &[u8]) -> Vec<u8> {
    let mut value = value.split(|&b| b == 0).next().unwrap_or_default().to_vec();
    value.push(0);
    value
}
fn c_atoll_bytes(value: &[u8]) -> i64 {
    unsafe { libc::atoll(c_input(value).as_ptr().cast()) as i64 }
}
fn c_atof_bytes(value: &[u8]) -> f64 {
    unsafe { libc::atof(c_input(value).as_ptr().cast()) }
}
macro_rules! impl_numeric {
    ($t:ty,$parse:expr) => {
        impl ReadOption for $t {
            fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
                Self::read_option_bytes(
                    dst,
                    &v.iter().map(|s| s.as_bytes().to_vec()).collect::<Vec<_>>(),
                )
                .map_err(|e| e.to_string())
            }
            fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
                *dst = ($parse)(&v[0]) as $t;
                Ok(())
            }
        }
    };
}
impl_numeric!(i32, |v: &[u8]| unsafe {
    libc::atoi(c_input(v).as_ptr().cast())
});
impl_numeric!(u32, |v: &[u8]| unsafe {
    libc::atoi(c_input(v).as_ptr().cast())
});
impl_numeric!(i64, c_atoll_bytes);
impl_numeric!(u64, c_atoll_bytes);
impl_numeric!(usize, c_atoll_bytes);
impl_numeric!(isize, c_atoll_bytes);
impl_numeric!(f64, c_atof_bytes);

impl ReadOption for bool {
    fn read_option(dst: &mut Self, _: &[String]) -> Result<(), String> {
        *dst = true;
        Ok(())
    }
    fn read_option_bytes(dst: &mut Self, _: &[Vec<u8>]) -> std::io::Result<()> {
        *dst = true;
        Ok(())
    }
    fn check_pcount(v: &[String], _: i32) -> bool {
        v.is_empty()
    }
    fn check_pcount_len(n: usize, _: i32) -> bool {
        n == 0
    }
}
impl ReadOption for String {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        *dst = v[0].clone();
        Ok(())
    }
    fn check_present(&self) -> bool {
        !self.is_empty()
    }
}
impl ReadOption for Vec<String> {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        *dst = v.to_vec();
        Ok(())
    }
    fn check_pcount(v: &[String], m: i32) -> bool {
        v.len() >= m as usize
    }
    fn check_pcount_len(n: usize, m: i32) -> bool {
        n >= m as usize
    }
}
impl ReadOption for Vec<u8> {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        *dst = v[0].as_bytes().to_vec();
        Ok(())
    }
    fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
        *dst = v[0].clone();
        Ok(())
    }
    fn check_present(&self) -> bool {
        !self.is_empty()
    }
}
impl ReadOption for Vec<Vec<u8>> {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        *dst = v.iter().map(|s| s.as_bytes().to_vec()).collect();
        Ok(())
    }
    fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
        *dst = v.to_vec();
        Ok(())
    }
    fn check_pcount(v: &[String], m: i32) -> bool {
        v.len() >= m as usize
    }
    fn check_pcount_len(n: usize, m: i32) -> bool {
        n >= m as usize
    }
}
impl ReadOption for std::ffi::OsString {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        *dst = v[0].clone().into();
        Ok(())
    }
    fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
        *dst = source_bytes_os(&v[0])?;
        Ok(())
    }
    fn read_option_os(dst: &mut Self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
        *dst = v[0].clone();
        Ok(())
    }
    fn check_present(&self) -> bool {
        !self.is_empty()
    }
}
impl ReadOption for std::path::PathBuf {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        *dst = v[0].clone().into();
        Ok(())
    }
    fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
        *dst = source_bytes_os(&v[0])?.into();
        Ok(())
    }
    fn read_option_os(dst: &mut Self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
        *dst = v[0].clone().into();
        Ok(())
    }
    fn check_present(&self) -> bool {
        !self.as_os_str().is_empty()
    }
}
macro_rules! impl_os_vector {
    ($t:ty) => {
        impl ReadOption for Vec<$t> {
            fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
                *dst = v.iter().cloned().map(Into::into).collect();
                Ok(())
            }
            fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
                *dst = v
                    .iter()
                    .map(|s| source_bytes_os(s).map(Into::into))
                    .collect::<std::io::Result<_>>()?;
                Ok(())
            }
            fn read_option_os(dst: &mut Self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
                *dst = v.iter().cloned().map(Into::into).collect();
                Ok(())
            }
            fn check_pcount(v: &[String], m: i32) -> bool {
                v.len() >= m as usize
            }
            fn check_pcount_len(n: usize, m: i32) -> bool {
                n >= m as usize
            }
        }
    };
}
impl_os_vector!(std::ffi::OsString);
impl_os_vector!(std::path::PathBuf);
impl<T: ReadOption + Default> ReadOption for OptionValue<T> {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        let mut value = T::default();
        T::read_option(&mut value, v)?;
        dst.set(value);
        Ok(())
    }
    fn read_option_bytes(dst: &mut Self, v: &[Vec<u8>]) -> std::io::Result<()> {
        let mut value = T::default();
        T::read_option_bytes(&mut value, v)?;
        dst.set(value);
        Ok(())
    }
    fn read_option_os(dst: &mut Self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
        let mut value = T::default();
        T::read_option_os(&mut value, v)?;
        dst.set(value);
        Ok(())
    }
    fn check_pcount(v: &[String], m: i32) -> bool {
        T::check_pcount(v, m)
    }
    fn check_pcount_len(n: usize, m: i32) -> bool {
        T::check_pcount_len(n, m)
    }
    fn check_present(&self) -> bool {
        self.present()
    }
    fn set_option_default(_: &mut Self, _: &Self) {}
    fn set_base_ptr(&mut self, base: &OptionBase) {
        let value = self.get_present().ok();
        *self = OptionValue::with_base(crate::util::options::OptionBase::new(
            &base.id,
            base.short_id,
            &base.desc,
            base.disabled,
            Some(&base.group_title),
        ));
        if let Some(value) = value {
            self.set(value);
        }
    }
}

pub fn read_option<T: ReadOption>(dst: &mut T, v: &[String]) -> Result<(), String> {
    T::read_option(dst, v)
}

pub fn read_option_bytes<T: ReadOption>(dst: &mut T, v: &[Vec<u8>]) -> std::io::Result<()> {
    T::read_option_bytes(dst, v)
}
pub fn read_option_os<T: ReadOption>(dst: &mut T, v: &[std::ffi::OsString]) -> std::io::Result<()> {
    T::read_option_os(dst, v)
}

pub fn check_pcount<T: ReadOption>(v: &[String], min_count: i32) -> bool {
    T::check_pcount(v, min_count)
}

pub fn check_present<T: ReadOption>(v: &T) -> bool {
    v.check_present()
}

pub fn set_option_default<T: ReadOption>(option: &mut T, value: &T) {
    T::set_option_default(option, value);
}

pub fn set_base_ptr<T: ReadOption>(option: &mut T, base: &OptionBase) {
    option.set_base_ptr(base);
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OptionBase {
    pub id: String,
    pub short_id: char,
    pub desc: String,
    pub disabled: bool,
    pub commands: Vec<u32>,
    pub group_title: String,
}

pub trait OptionDescBase {
    fn base(&self) -> &OptionBase;
    fn read(&self, v: &[String]) -> Result<(), String>;
    fn read_bytes(&self, v: &[Vec<u8>]) -> std::io::Result<()> {
        let values = v
            .iter()
            .map(|s| String::from_utf8(s.clone()).map_err(std::io::Error::other))
            .collect::<std::io::Result<Vec<_>>>()?;
        self.read(&values)
            .map_err(|message| raw_error(message.into_bytes()))
    }
    fn read_os(&self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
        let values = v
            .iter()
            .map(|s| source_os_bytes(s).map(|s| s.into_owned()))
            .collect::<std::io::Result<Vec<_>>>()?;
        self.read_bytes(&values)
    }
    fn present(&self) -> bool;
    fn set_default(&self);
}

pub struct OptionDesc<T: ReadOption> {
    base: OptionBase,
    default: T,
    min_count: i32,
    store: Rc<RefCell<T>>,
}

impl<T: ReadOption> OptionDesc<T> {
    pub fn new(
        id: &str,
        short_id: char,
        desc: &str,
        disabled: bool,
        store: Rc<RefCell<T>>,
        default: T,
        min_count: i32,
        group_title: &str,
        commands: &[u32],
    ) -> Self {
        let base = OptionBase {
            id: id.to_string(),
            short_id,
            desc: desc.to_string(),
            disabled,
            commands: commands.to_vec(),
            group_title: group_title.to_string(),
        };
        store.borrow_mut().set_base_ptr(&base);
        Self {
            base,
            default,
            min_count,
            store,
        }
    }
}

impl<T: ReadOption> OptionDesc<T> {
    fn validate_count(&self, n: usize) -> std::io::Result<()> {
        if T::check_pcount_len(n, self.min_count) {
            return Ok(());
        }
        let prefix = if self.base.short_id != '\0' {
            format!("-{}/", self.base.short_id)
        } else {
            String::new()
        };
        Err(raw_error(
            format!(
                "Invalid parameter count for option '{}--{}'",
                prefix, self.base.id
            )
            .into_bytes(),
        ))
    }
}

impl<T: ReadOption> OptionDescBase for OptionDesc<T> {
    fn base(&self) -> &OptionBase {
        &self.base
    }

    fn read(&self, v: &[String]) -> Result<(), String> {
        self.read_bytes(&v.iter().map(|s| s.as_bytes().to_vec()).collect::<Vec<_>>())
            .map_err(|e| e.to_string())
    }

    fn read_bytes(&self, v: &[Vec<u8>]) -> std::io::Result<()> {
        self.validate_count(v.len())?;
        T::read_option_bytes(&mut self.store.borrow_mut(), v)
    }
    fn read_os(&self, v: &[std::ffi::OsString]) -> std::io::Result<()> {
        self.validate_count(v.len())?;
        T::read_option_os(&mut self.store.borrow_mut(), v)
    }
    fn present(&self) -> bool {
        self.store.borrow().check_present()
    }

    fn set_default(&self) {
        T::set_option_default(&mut self.store.borrow_mut(), &self.default);
    }
}

pub type OptionRef = Rc<dyn OptionDescBase>;

#[derive(Clone)]
pub struct OptionsGroup {
    pub options: Vec<OptionRef>,
    pub title: String,
    pub commands: Vec<u32>,
    pub disabled: bool,
}

impl OptionsGroup {
    pub fn new(title: &str, commands: Vec<u32>, disabled: bool) -> Self {
        Self {
            options: Vec::new(),
            title: title.to_string(),
            commands,
            disabled,
        }
    }
}

pub struct CommandLineParser {
    map: BTreeMap<String, OptionRef>,
    map_short: BTreeMap<char, OptionRef>,
    groups: Vec<OptionsGroup>,
    command_codes: BTreeMap<String, u32>,
    commands: Vec<(String, String)>,
}

impl Default for CommandLineParser {
    fn default() -> Self {
        Self::new()
    }
}

impl CommandLineParser {
    pub fn new() -> Self {
        Self {
            map: BTreeMap::new(),
            map_short: BTreeMap::new(),
            groups: Vec::new(),
            command_codes: BTreeMap::new(),
            commands: Vec::new(),
        }
    }

    pub fn add_group(&mut self, title: &str, commands: Vec<u32>, disabled: bool) -> usize {
        self.groups
            .push(OptionsGroup::new(title, commands, disabled));
        self.groups.len() - 1
    }

    pub fn add_option<T: ReadOption>(
        &mut self,
        group_idx: usize,
        id: &str,
        short_id: char,
        desc: &str,
        store: Rc<RefCell<T>>,
        default: T,
        min_count: i32,
    ) -> &mut Self {
        let group = &self.groups[group_idx];
        let option: OptionRef = Rc::new(OptionDesc::new(
            id,
            short_id,
            desc,
            group.disabled,
            store,
            default,
            min_count,
            &group.title,
            &group.commands,
        ));
        self.groups[group_idx].options.push(option.clone());
        self.map.insert(id.to_string(), option.clone());
        self.map_short.insert(short_id, option);
        self
    }

    pub fn add_command(&mut self, s: &str, desc: &str, code: u32) -> &mut Self {
        self.command_codes.insert(s.to_string(), code);
        self.commands.push((s.to_string(), desc.to_string()));
        self
    }

    pub fn store(&self, args: &[&str], command: &mut u32) -> Result<(), String> {
        self.store_bytes(
            &args.iter().map(|s| s.as_bytes()).collect::<Vec<_>>(),
            command,
        )
        .map_err(|e| e.to_string())
    }
    pub fn store_bytes(&self, args: &[&[u8]], command: &mut u32) -> std::io::Result<()> {
        self.store_tokens(
            args.iter().map(|s| Argument::Bytes(s.to_vec())).collect(),
            command,
        )
    }
    pub fn store_os(&self, args: &[&std::ffi::OsStr], command: &mut u32) -> std::io::Result<()> {
        self.store_tokens(
            args.iter().map(|s| Argument::Os((*s).to_owned())).collect(),
            command,
        )
    }
    fn store_tokens(&self, args: Vec<Argument>, command: &mut u32) -> std::io::Result<()> {
        if args.len() < 2 {
            return Err(raw_error(
                b"Syntax: diamond COMMAND [OPTIONS]. To print help message: diamond help".to_vec(),
            ));
        }
        let cmd = args[1].bytes()?;
        let cmd = cmd.strip_prefix(b"--").unwrap_or(&cmd);
        *command = *std::str::from_utf8(cmd)
            .ok()
            .and_then(|s| self.command_codes.get(s))
            .ok_or_else(|| {
                let mut message = b"Invalid command: ".to_vec();
                message.extend_from_slice(cmd);
                message.extend_from_slice(b". To print help message: diamond help");
                raw_error(message)
            })?;
        for option in self.map.values() {
            option.set_default();
        }
        let mut group = Vec::new();
        for arg in args.into_iter().skip(2) {
            let bytes = arg.bytes()?;
            if bytes.len() > 1 && bytes[0] == b'-' && !bytes[1].is_ascii_digit() {
                self.store_token_option(&group, *command)?;
                group.clear();
            }
            group.push(arg);
        }
        self.store_token_option(&group, *command)
    }
    fn store_token_option(&self, v: &[Argument], command: u32) -> std::io::Result<()> {
        if v.is_empty() {
            return Ok(());
        }
        let token = v[0].bytes()?;
        if token.len() <= 1 {
            return Err(raw_error(b"Invalid option syntax.".to_vec()));
        }
        let (id, option, attached) = if token.starts_with(b"--") {
            let id = &token[2..];
            (
                id,
                std::str::from_utf8(id)
                    .ok()
                    .and_then(|s| self.map.get(s))
                    .cloned(),
                false,
            )
        } else if token[0] == b'-' {
            (
                &token[1..2],
                self.map_short.get(&(token[1] as char)).cloned(),
                token.len() > 2,
            )
        } else {
            return Err(raw_error(
                b"Command line options must begin with - or --.".to_vec(),
            ));
        };
        let error = |prefix: &[u8]| {
            let mut message = prefix.to_vec();
            message.extend_from_slice(id);
            raw_error(message)
        };
        let option = option.ok_or_else(|| error(b"Invalid option: "))?;
        if option.base().disabled {
            return Err(error(b"Invalid option: "));
        }
        if !cfg!(feature = "extra") && !option.base().commands.contains(&command) {
            return Err(error(b"Option is not permitted for this workflow: "));
        }
        let mut values = Vec::new();
        if attached {
            values.push(v[0].attached()?);
        }
        values.extend_from_slice(&v[1..]);
        if matches!(v[0], Argument::Os(_)) {
            let values = values
                .into_iter()
                .map(|v| match v {
                    Argument::Os(s) => s,
                    Argument::Bytes(_) => unreachable!(),
                })
                .collect::<Vec<_>>();
            option.read_os(&values)
        } else {
            let values = values
                .into_iter()
                .map(|v| match v {
                    Argument::Bytes(s) => s,
                    Argument::Os(_) => unreachable!(),
                })
                .collect::<Vec<_>>();
            option.read_bytes(&values)
        }
    }

    pub fn require(&self, option: &str) -> Result<(), String> {
        let o = self
            .map
            .get(option)
            .ok_or_else(|| "Unknown option.".to_string())?;
        if !o.present() {
            return Err(format!(
                "Missing parameter: --{}/-{}",
                o.base().id,
                o.base().short_id.to_string().trim_end_matches('\0')
            ));
        }
        Ok(())
    }

    pub fn print_help(&self) -> String {
        const COL1_WIDTH: usize = 25;
        let mut out = String::new();
        out.push_str("Syntax: diamond COMMAND [OPTIONS]\n\nCommands:\n");
        for (cmd, desc) in &self.commands {
            if desc.is_empty() {
                continue;
            }
            out.push_str(cmd);
            for _ in 0..COL1_WIDTH.saturating_sub(cmd.len()) {
                out.push(' ');
            }
            out.push_str(desc);
            out.push('\n');
        }
        out.push_str(
            "\nPossible [OPTIONS] for COMMAND can be seen with syntax: diamond COMMAND\n\n",
        );
        out.push_str("Online documentation at http://www.diamondsearch.org\n");
        out
    }

    pub fn print_documentation(&self, command: u32) -> String {
        const COL1_WIDTH: usize = 25;
        let mut out = String::from("Options:\n");
        for group in &self.groups {
            // The C++ implementation intentionally iterates the command list,
            // so a duplicated command entry prints the group more than once.
            for &group_command in &group.commands {
                if group_command != command {
                    continue;
                }
                for option in &group.options {
                    if option.base().desc.is_empty() {
                        continue;
                    }
                    let mut col1 = format!("--{}", option.base().id);
                    let pad = COL1_WIDTH.max(col1.len()) - col1.len();
                    col1.extend(std::iter::repeat_n(' ', pad));
                    out.push_str(&col1);
                    out.push_str(&option.base().desc);
                    out.push('\n');
                }
            }
        }
        out.push('\n');
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_values_numeric_crt_and_optionvalue_share_parser() {
        let text = Rc::new(RefCell::new(Vec::<u8>::new()));
        let values = Rc::new(RefCell::new(OptionValue::<Vec<Vec<u8>>>::new()));
        let floating = Rc::new(RefCell::new(0.0f64));
        let mut parser = CommandLineParser::new();
        parser.add_command("run", "", 1);
        let group = parser.add_group("test", vec![1], false);
        parser
            .add_option(group, "text", 't', "", text.clone(), Vec::<u8>::new(), 1)
            .add_option(
                group,
                "values",
                'v',
                "",
                values.clone(),
                OptionValue::new(),
                1,
            )
            .add_option(group, "float", 'f', "", floating.clone(), 0.0, 1);
        let mut command = 0;
        parser
            .store_bytes(
                &[
                    b"diamond",
                    b"run",
                    b"-t\xff",
                    b"--values",
                    b"\xfe",
                    b"a",
                    b"--values",
                    b"\xfd",
                    b"--float",
                    b"0x1.8p+2tail",
                ],
                &mut command,
            )
            .unwrap();
        assert_eq!(*text.borrow(), b"\xff");
        assert_eq!(
            values.borrow().get_present().unwrap(),
            vec![b"\xfd".to_vec()]
        );
        assert_eq!(*floating.borrow(), 6.0);
        let mut integer = 99i64;
        ReadOption::read_option_bytes(&mut integer, &[b"\xc2\xa012".to_vec()]).unwrap();
        assert_eq!(integer, 0); // Unicode NBSP is not a C-locale space.
        let error = parser
            .store_bytes(&[b"diamond", b"run", b"--\xff"], &mut command)
            .unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(b"Invalid option: \xff".as_slice())
        );
    }

    #[cfg(unix)]
    #[test]
    fn os_paths_keep_raw_bytes_in_scalar_vector_and_optional_cells() {
        use std::ffi::{OsStr, OsString};
        use std::os::unix::ffi::OsStringExt;
        use std::path::PathBuf;
        let path = Rc::new(RefCell::new(PathBuf::new()));
        let paths = Rc::new(RefCell::new(OptionValue::<Vec<PathBuf>>::new()));
        let mut parser = CommandLineParser::new();
        parser.add_command("run", "", 1);
        let group = parser.add_group("test", vec![1], false);
        parser
            .add_option(group, "path", 'p', "", path.clone(), PathBuf::new(), 1)
            .add_option(
                group,
                "paths",
                'q',
                "",
                paths.clone(),
                OptionValue::new(),
                1,
            );
        let attached = OsString::from_vec(b"-pdir/\xff".to_vec());
        let second = OsString::from_vec(b"other/\xfe".to_vec());
        parser
            .store_os(
                &[
                    OsStr::new("diamond"),
                    OsStr::new("run"),
                    attached.as_os_str(),
                    OsStr::new("--paths"),
                    second.as_os_str(),
                ],
                &mut 0,
            )
            .unwrap();
        assert_eq!(
            *path.borrow(),
            PathBuf::from(OsString::from_vec(b"dir/\xff".to_vec()))
        );
        assert_eq!(
            paths.borrow().get_present().unwrap(),
            vec![PathBuf::from(second)]
        );
    }

    #[cfg(windows)]
    #[test]
    fn os_path_cells_retain_utf16_even_when_source_narrow_text_cannot() {
        use std::ffi::{OsStr, OsString};
        use std::os::windows::ffi::{OsStrExt, OsStringExt};
        use std::path::PathBuf;
        let path = Rc::new(RefCell::new(PathBuf::new()));
        let mut parser = CommandLineParser::new();
        parser.add_command("run", "", 1);
        let group = parser.add_group("test", vec![1], false);
        parser.add_option(group, "path", 'p', "", path.clone(), PathBuf::new(), 1);
        let raw = OsString::from_wide(&[b'x' as u16, 0xd800, b'y' as u16]);
        parser
            .store_os(
                &[
                    OsStr::new("diamond"),
                    OsStr::new("run"),
                    OsStr::new("--path"),
                    raw.as_os_str(),
                ],
                &mut 0,
            )
            .unwrap();
        assert_eq!(
            path.borrow().as_os_str().encode_wide().collect::<Vec<_>>(),
            vec![b'x' as u16, 0xd800, b'y' as u16]
        );
    }

    #[test]
    fn test_read_option_helpers() {
        let mut x = 0i32;
        read_option(&mut x, &["7".to_string()]).unwrap();
        assert_eq!(x, 7);
        assert!(check_pcount::<i32>(&["x".to_string()], 1));
        assert!(check_pcount::<bool>(&[], 1));
        assert!(check_pcount::<Vec<String>>(
            &["a".to_string(), "b".to_string()],
            2
        ));

        let mut s = String::new();
        assert!(!check_present(&s));
        read_option(&mut s, &["db".to_string()]).unwrap();
        assert!(check_present(&s));
        set_option_default(&mut s, &"x".to_string());
        assert_eq!(s, "x");

        let mut opt_vec = OptionValue::<Vec<String>>::new();
        read_option(&mut opt_vec, &["a".to_string(), "b".to_string()]).unwrap();
        assert_eq!(
            opt_vec.get_present().unwrap(),
            vec!["a".to_string(), "b".to_string()]
        );
        assert!(check_pcount::<OptionValue<Vec<String>>>(
            &["a".to_string(), "b".to_string()],
            2
        ));

        let mut signed = 0i32;
        read_option(&mut signed, &["  -12trailing".to_string()]).unwrap();
        assert_eq!(signed, -12);
        let mut unsigned = 0u32;
        read_option(&mut unsigned, &["-1".to_string()]).unwrap();
        assert_eq!(unsigned, u32::MAX);
        let mut floating = 0.0f64;
        read_option(&mut floating, &["2.5e1suffix".to_string()]).unwrap();
        assert_eq!(floating, 25.0);
        let mut optional = OptionValue::<i64>::new();
        read_option(&mut optional, &["not-a-number".to_string()]).unwrap();
        assert_eq!(optional.get_present(), Ok(0));
    }

    #[test]
    fn test_command_line_parser_store_require_and_docs() {
        let db = Rc::new(RefCell::new(String::new()));
        let threads = Rc::new(RefCell::new(1u32));
        let verbose = Rc::new(RefCell::new(false));
        let mut parser = CommandLineParser::new();
        parser.add_command("blastp", "protein search", 1);
        parser.add_command("makedb", "", 2);
        let group = parser.add_group("General", vec![1], false);
        parser
            .add_option(group, "db", 'd', "database", db.clone(), String::new(), 1)
            .add_option(group, "threads", 'p', "threads", threads.clone(), 1u32, 1)
            .add_option(group, "verbose", 'v', "verbose", verbose.clone(), false, 1);

        let mut command = 0u32;
        parser
            .store(
                &["diamond", "blastp", "--db", "nr", "-p8", "--verbose"],
                &mut command,
            )
            .unwrap();
        assert_eq!(command, 1);
        assert_eq!(db.borrow().as_str(), "nr");
        assert_eq!(*threads.borrow(), 8);
        assert!(*verbose.borrow());
        parser.require("db").unwrap();

        let help = parser.print_help();
        assert!(help.contains("blastp"));
        assert!(!help.contains("makedb                         "));
        let docs = parser.print_documentation(1);
        assert!(docs.contains("--threads"));

        let err = parser.store(&["diamond", "makedb", "--db", "x"], &mut command);
        if cfg!(feature = "extra") {
            err.unwrap();
            assert_eq!(db.borrow().as_str(), "x");
        } else {
            assert_eq!(
                err.unwrap_err(),
                "Option is not permitted for this workflow: db"
            );
            assert!(db.borrow().is_empty());
        }
        // store resets every option before applying the next command, even
        // when EXTRA bypasses group membership or baseline rejects the store.
        assert_eq!(command, 2);
        assert_eq!(*threads.borrow(), 1);
        assert!(!*verbose.borrow());
    }

    #[test]
    fn test_command_line_parser_exact_help_and_documentation() {
        let value = Rc::new(RefCell::new(String::new()));
        let mut parser = CommandLineParser::new();
        parser
            .add_command("blastp", "protein search", 1)
            .add_command("hidden", "", 2);
        let group = parser.add_group("General", vec![1, 1], false);
        parser.add_option(
            group,
            "database",
            'd',
            "database path",
            value,
            String::new(),
            1,
        );

        assert_eq!(
            parser.print_help(),
            "Syntax: diamond COMMAND [OPTIONS]\n\nCommands:\nblastp                   protein search\n\nPossible [OPTIONS] for COMMAND can be seen with syntax: diamond COMMAND\n\nOnline documentation at http://www.diamondsearch.org\n"
        );
        assert_eq!(
            parser.print_documentation(1),
            "Options:\n--database               database path\n--database               database path\n\n"
        );
    }

    #[test]
    fn test_command_line_parser_exact_error_paths_and_default_reset() {
        let name = Rc::new(RefCell::new(String::new()));
        let disabled = Rc::new(RefCell::new(false));
        let optional = Rc::new(RefCell::new(OptionValue::<String>::new()));
        let mut parser = CommandLineParser::new();
        parser.add_command("run", "run", 7);
        let enabled_group = parser.add_group("Enabled", vec![7], false);
        parser.add_option(
            enabled_group,
            "name",
            'n',
            "name",
            name.clone(),
            "default".to_string(),
            1,
        );
        parser.add_option(
            enabled_group,
            "optional",
            'o',
            "optional",
            optional.clone(),
            OptionValue::new(),
            1,
        );
        let disabled_group = parser.add_group("Disabled", vec![7], true);
        parser.add_option(
            disabled_group,
            "disabled",
            'x',
            "disabled",
            disabled,
            false,
            1,
        );

        let mut command = 0;
        assert_eq!(
            parser.store(&["diamond"], &mut command).unwrap_err(),
            "Syntax: diamond COMMAND [OPTIONS]. To print help message: diamond help"
        );
        assert_eq!(
            parser.store(&["diamond", "bad"], &mut command).unwrap_err(),
            "Invalid command: bad. To print help message: diamond help"
        );
        assert_eq!(
            parser
                .store(&["diamond", "run", "value"], &mut command)
                .unwrap_err(),
            "Command line options must begin with - or --."
        );
        assert_eq!(
            parser
                .store(&["diamond", "run", "-"], &mut command)
                .unwrap_err(),
            "Invalid option syntax."
        );
        assert_eq!(
            parser
                .store(&["diamond", "run", "--missing"], &mut command)
                .unwrap_err(),
            "Invalid option: missing"
        );
        assert_eq!(
            parser
                .store(&["diamond", "run", "--disabled"], &mut command)
                .unwrap_err(),
            "Invalid option: disabled"
        );
        assert_eq!(
            parser
                .store(&["diamond", "run", "--name"], &mut command)
                .unwrap_err(),
            "Invalid parameter count for option '-n/--name'"
        );
        assert_eq!(parser.require("unknown").unwrap_err(), "Unknown option.");
        assert_eq!(
            optional.borrow().require().unwrap_err(),
            "Missing parameter: --optional/-o"
        );

        parser
            .store(&["diamond", "--run", "-ncustom"], &mut command)
            .unwrap();
        assert_eq!(command, 7);
        assert_eq!(name.borrow().as_str(), "custom");
        parser.store(&["diamond", "run"], &mut command).unwrap();
        assert_eq!(name.borrow().as_str(), "default");
    }
}
