use crate::basic::value::TaxId;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind};

pub fn read_nodes_dmp<F>(file_name: &str, mut f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    let file = File::open(file_name)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            break;
        }
        let mut fields = line.split("\t|\t");
        let taxid = fields
            .next()
            .ok_or_else(|| Error::new(ErrorKind::InvalidData, "missing taxid"))?
            .trim()
            .parse::<TaxId>()
            .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
        let parent = fields
            .next()
            .ok_or_else(|| Error::new(ErrorKind::InvalidData, "missing parent taxid"))?
            .trim()
            .parse::<TaxId>()
            .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
        let rank = fields.next().unwrap_or("");
        f(taxid, parent, rank);
    }
    Ok(())
}

pub fn read_names_dmp<F>(file_name: &str, mut f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    let file = File::open(file_name)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split("\t|\t");
        let id = fields
            .next()
            .ok_or_else(|| Error::new(ErrorKind::InvalidData, "missing taxid"))?
            .trim()
            .parse::<TaxId>()
            .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
        let name = fields.next().unwrap_or("");
        let _unique_name = fields.next();
        let type_ = fields.next().unwrap_or("").trim_end_matches(['\t', '|']);
        if type_ == "scientific name" {
            f(id, name);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_nodes_dmp_template_callback() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-nodes-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "1\t|\t1\t|\tno rank\t|\t").unwrap();
            writeln!(f, "2\t|\t1\t|\tsuperkingdom\t|\t").unwrap();
        }

        let mut rows = Vec::new();
        read_nodes_dmp(path.to_str().unwrap(), |taxid, parent, rank| {
            rows.push((taxid, parent, rank.to_string()));
        })
        .unwrap();
        assert_eq!(
            rows,
            vec![
                (1, 1, "no rank".to_string()),
                (2, 1, "superkingdom".to_string())
            ]
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_nodes_dmp_stops_at_empty_line() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-nodes-empty-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "1\t|\t1\t|\tno rank\t|\t").unwrap();
            writeln!(f).unwrap();
            writeln!(f, "2\t|\t1\t|\tsuperkingdom\t|\t").unwrap();
        }

        let mut rows = Vec::new();
        read_nodes_dmp(path.to_str().unwrap(), |taxid, parent, rank| {
            rows.push((taxid, parent, rank.to_string()));
        })
        .unwrap();
        assert_eq!(rows, vec![(1, 1, "no rank".to_string())]);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_nodes_dmp_rejects_malformed_taxid() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-nodes-bad-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "not-a-taxid\t|\t1\t|\tno rank\t|\t").unwrap();
        }

        let err = read_nodes_dmp(path.to_str().unwrap(), |_, _, _| {}).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_names_dmp_filters_scientific_name() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-names-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "2\t|\tBacteria\t|\t\t|\tscientific name\t|").unwrap();
            writeln!(f, "2\t|\tbacteria\t|\t\t|\tcommon name\t|").unwrap();
            writeln!(f, "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|").unwrap();
        }

        let mut rows = Vec::new();
        read_names_dmp(path.to_str().unwrap(), |id, name| {
            rows.push((id, name.to_string()));
        })
        .unwrap();
        assert_eq!(
            rows,
            vec![
                (2, "Bacteria".to_string()),
                (9606, "Homo sapiens".to_string())
            ]
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_names_dmp_rejects_malformed_taxid() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-names-bad-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "bad\t|\tBacteria\t|\t\t|\tscientific name\t|").unwrap();
        }

        let err = read_names_dmp(path.to_str().unwrap(), |_, _| {}).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        std::fs::remove_file(path).unwrap();
    }
}
