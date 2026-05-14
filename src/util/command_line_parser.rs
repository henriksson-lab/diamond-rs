use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::util::options::OptionValue;

pub trait ReadOption: Clone + 'static {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String>;
    fn check_pcount(v: &[String], min_count: i32) -> bool {
        let _ = min_count;
        v.len() == 1
    }
    fn check_present(&self) -> bool {
        false
    }
    fn set_option_default(dst: &mut Self, value: &Self) {
        *dst = value.clone();
    }
}

macro_rules! impl_read_option_from_str {
    ($($t:ty),* $(,)?) => {
        $(
            impl ReadOption for $t {
                fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
                    *dst = v[0].parse::<$t>().unwrap_or_default();
                    Ok(())
                }
            }
        )*
    };
}

impl_read_option_from_str!(i32, u32, i64, u64, usize, isize, f64);

impl ReadOption for bool {
    fn read_option(dst: &mut Self, _v: &[String]) -> Result<(), String> {
        *dst = true;
        Ok(())
    }

    fn check_pcount(v: &[String], _min_count: i32) -> bool {
        v.is_empty()
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

    fn check_pcount(v: &[String], min_count: i32) -> bool {
        v.len() >= min_count as usize
    }
}

macro_rules! impl_read_option_value_from_str {
    ($($t:ty),* $(,)?) => {
        $(
            impl ReadOption for OptionValue<$t> {
                fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
                    let value = v[0]
                        .parse::<$t>()
                        .map_err(|_| "Invalid option value.".to_string())?;
                    dst.set(value);
                    Ok(())
                }

                fn check_present(&self) -> bool {
                    self.present()
                }

                fn set_option_default(_dst: &mut Self, _value: &Self) {}
            }
        )*
    };
}

impl_read_option_value_from_str!(i32, u32, i64, u64, usize, isize, f64, String);

impl ReadOption for OptionValue<Vec<String>> {
    fn read_option(dst: &mut Self, v: &[String]) -> Result<(), String> {
        dst.set(v.to_vec());
        Ok(())
    }

    fn check_pcount(v: &[String], min_count: i32) -> bool {
        v.len() >= min_count as usize
    }

    fn check_present(&self) -> bool {
        self.present()
    }

    fn set_option_default(_dst: &mut Self, _value: &Self) {}
}

pub fn read_option<T: ReadOption>(dst: &mut T, v: &[String]) -> Result<(), String> {
    T::read_option(dst, v)
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
        Self {
            base: OptionBase {
                id: id.to_string(),
                short_id,
                desc: desc.to_string(),
                disabled,
                commands: commands.to_vec(),
                group_title: group_title.to_string(),
            },
            default,
            min_count,
            store,
        }
    }
}

impl<T: ReadOption> OptionDescBase for OptionDesc<T> {
    fn base(&self) -> &OptionBase {
        &self.base
    }

    fn read(&self, v: &[String]) -> Result<(), String> {
        if !T::check_pcount(v, self.min_count) {
            let prefix = if self.base.short_id != '\0' {
                format!("-{}/", self.base.short_id)
            } else {
                String::new()
            };
            return Err(format!(
                "Invalid parameter count for option '{}--{}'",
                prefix, self.base.id
            ));
        }
        T::read_option(&mut self.store.borrow_mut(), v)
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
        if args.len() < 2 {
            return Err(
                "Syntax: diamond COMMAND [OPTIONS]. To print help message: diamond help"
                    .to_string(),
            );
        }
        let mut cmd = args[1].to_string();
        if cmd.starts_with("--") {
            cmd = cmd[2..].to_string();
        }
        *command = *self.command_codes.get(&cmd).ok_or_else(|| {
            format!("Invalid command: {cmd}. To print help message: diamond help")
        })?;

        for option in self.map.values() {
            option.set_default();
        }

        let mut v = Vec::new();
        for arg in args.iter().skip(2) {
            if arg.starts_with('-') && arg.len() > 1 && !arg.as_bytes()[1].is_ascii_digit() {
                self.store_option(&v, *command)?;
                v.clear();
            }
            v.push((*arg).to_string());
        }
        self.store_option(&v, *command)
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
                o.base().short_id
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
            if !group.commands.contains(&command) {
                continue;
            }
            for option in &group.options {
                if option.base().desc.is_empty() {
                    continue;
                }
                let mut col1 = format!("--{}", option.base().id);
                let pad = COL1_WIDTH.max(col1.len()) - col1.len();
                col1.extend(std::iter::repeat(' ').take(pad));
                out.push_str(&col1);
                out.push_str(&option.base().desc);
                out.push('\n');
            }
        }
        out.push('\n');
        out
    }

    fn store_option(&self, v: &[String], command: u32) -> Result<(), String> {
        if v.is_empty() {
            return Ok(());
        }
        let o;
        let id;
        let mut v2 = Vec::new();
        if v[0].len() <= 1 {
            return Err("Invalid option syntax.".to_string());
        } else if v[0].starts_with("--") {
            id = v[0][2..].to_string();
            o = self.map.get(&id).cloned();
        } else if v[0].starts_with('-') {
            id = v[0].chars().nth(1).unwrap().to_string();
            let short = v[0].chars().nth(1).unwrap();
            o = self.map_short.get(&short).cloned();
            if v[0].len() > 2 {
                v2.push(v[0][2..].to_string());
            }
        } else {
            return Err("Command line options must begin with - or --.".to_string());
        }

        let o = o.ok_or_else(|| format!("Invalid option: {id}"))?;
        if o.base().disabled {
            return Err(format!("Invalid option: {id}"));
        }
        if !o.base().commands.contains(&command) {
            return Err(format!("Option is not permitted for this workflow: {id}"));
        }
        v2.extend_from_slice(&v[1..]);
        o.read(&v2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(
            err.unwrap_err(),
            "Option is not permitted for this workflow: db"
        );
    }
}
