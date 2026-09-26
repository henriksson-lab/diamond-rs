//! BLAST diagnostic message ownership from
//! `diamond/src/lib/blast/blast_message.cpp`.

pub const BLAST_MESSAGE_NO_CONTEXT: i32 = -1;
pub const BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS: &str =
    "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum BlastSeverity {
    Info = 1,
    Warning = 2,
    Error = 3,
    Fatal = 4,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MessageOrigin {
    pub filename: String,
    pub line_number: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastMessage {
    pub next: Option<Box<BlastMessage>>,
    pub severity: BlastSeverity,
    pub message: String,
    pub origin: Option<MessageOrigin>,
    pub context: i32,
}

/// C++ `SMessageOriginNew`.
///
/// The vendored implementation initializes its local result pointer to null,
/// never allocates it, and therefore returns null for every input. This keeps
/// that active behavior rather than silently repairing an upstream bug.
pub fn s_message_origin_new(filename: Option<&str>, _line_number: u32) -> Option<MessageOrigin> {
    if filename.is_none_or(str::is_empty) {
        return None;
    }
    None
}

/// C++ `SMessageOriginFree`; ownership consumption is Rust's deallocation.
pub fn s_message_origin_free(_origin: Option<MessageOrigin>) -> Option<MessageOrigin> {
    None
}

/// C++ `Blast_MessageFree`; detach and drop every linked node and origin.
pub fn blast_message_free(mut message: Option<Box<BlastMessage>>) -> Option<Box<BlastMessage>> {
    while let Some(mut current) = message {
        message = current.next.take();
        let _ = s_message_origin_free(current.origin.take());
    }
    None
}

/// Append one message and return the C API status code.
///
/// A missing pointer-to-head returns `1`. Allocation failure while copying the
/// message returns `-1`; success returns `0` and preserves insertion order.
pub fn blast_message_write(
    head: Option<&mut Option<Box<BlastMessage>>>,
    severity: BlastSeverity,
    context: i32,
    message: &str,
) -> i32 {
    let Some(head) = head else {
        return 1;
    };
    let mut owned_message = String::new();
    if owned_message.try_reserve_exact(message.len()).is_err() {
        return -1;
    }
    owned_message.push_str(message);
    let new_message = Box::new(BlastMessage {
        next: None,
        severity,
        message: owned_message,
        origin: None,
        context,
    });
    let mut cursor = head;
    while let Some(message) = cursor {
        cursor = &mut message.next;
    }
    *cursor = Some(new_message);
    0
}

/// C++ `Blast_MessagePost`. The vendored body validates only the pointer and
/// performs no logging or traversal.
pub const fn blast_message_post(message: Option<&BlastMessage>) -> i32 {
    if message.is_some() {
        0
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_and_severity_values_match_the_c_api() {
        assert_eq!(BLAST_MESSAGE_NO_CONTEXT, -1);
        assert_eq!(BlastSeverity::Info as i32, 1);
        assert_eq!(BlastSeverity::Warning as i32, 2);
        assert_eq!(BlastSeverity::Error as i32, 3);
        assert_eq!(BlastSeverity::Fatal as i32, 4);
        assert_eq!(
            BLAST_ERR_MSG_CANT_CALCULATE_UNGAPPED_KA_PARAMS,
            "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options"
        );
    }

    #[test]
    fn origin_constructor_preserves_the_vendored_null_result_bug() {
        assert_eq!(s_message_origin_new(None, 4), None);
        assert_eq!(s_message_origin_new(Some(""), 4), None);
        assert_eq!(s_message_origin_new(Some("blast.c"), 4), None);
        assert_eq!(
            s_message_origin_free(Some(MessageOrigin {
                filename: "blast.c".into(),
                line_number: 4,
            })),
            None
        );
    }

    #[test]
    fn write_rejects_a_missing_head_and_appends_in_order() {
        assert_eq!(
            blast_message_write(None, BlastSeverity::Error, 1, "ignored"),
            1
        );
        let mut head = None;
        assert_eq!(
            blast_message_write(
                Some(&mut head),
                BlastSeverity::Warning,
                BLAST_MESSAGE_NO_CONTEXT,
                "first",
            ),
            0
        );
        assert_eq!(
            blast_message_write(Some(&mut head), BlastSeverity::Fatal, 7, "second"),
            0
        );
        let first = head.as_deref().unwrap();
        assert_eq!(first.severity, BlastSeverity::Warning);
        assert_eq!(first.context, BLAST_MESSAGE_NO_CONTEXT);
        assert_eq!(first.message, "first");
        assert_eq!(first.origin, None);
        let second = first.next.as_deref().unwrap();
        assert_eq!(second.severity, BlastSeverity::Fatal);
        assert_eq!(second.context, 7);
        assert_eq!(second.message, "second");
        assert!(second.next.is_none());
    }

    #[test]
    fn post_and_free_return_the_original_status_and_null_results() {
        assert_eq!(blast_message_post(None), 1);
        let message = Box::new(BlastMessage {
            next: None,
            severity: BlastSeverity::Info,
            message: "info".into(),
            origin: None,
            context: 0,
        });
        assert_eq!(blast_message_post(Some(&message)), 0);
        assert_eq!(blast_message_free(Some(message)), None);
        assert_eq!(blast_message_free(None), None);
    }
}
