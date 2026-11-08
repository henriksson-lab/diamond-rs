//! Gap-affine wavefront alignment configured by DIAMOND's DNA wrappers.
//!
//! Translated from DIAMOND's recovered WFA2-lib.diamond v2.3 tree
//! `3dc90a492f4a2397114e0c51bc168f1471dd398a`, last present at Source commit
//! `a9b7817292a075376858771031439f4dc913eae3`. Its DNA caller tree is identical
//! to pinned Source `1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`; later history
//! deletes the library without changing those callers. The active profile is
//! match cost zero, affine gaps, ASCII, high-memory traceback,
//! WFadaptive(10,50,1), XDrop(100,1). Canonical WFA2 was the initial porting
//! reference; recovered Source active-bound pruning/backtrace semantics apply.
//!
//! MIT License
//! Copyright (c) 2017 Santiago Marco-Sola
//!
//! Permission is hereby granted, free of charge, to any person obtaining a copy
//! of this software and associated documentation files (the "Software"), to deal
//! in the Software without restriction, including without limitation the rights
//! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! copies of the Software, and to permit persons to whom the Software is
//! furnished to do so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included in
//! all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//! SOFTWARE.

const NULL: i32 = i32::MIN / 2;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// WFA's completed (0) or partial/trimmed (1) status.
    pub status: i32,
    pub wavefront_score: i32,
    /// Internal matches are M; mismatches X; I consumes text, D pattern.
    pub operations: Vec<u8>,
}

impl Alignment {
    /// Source getCIGAR(true) writes matches as '=' and preserves mismatches.
    pub fn cigar_string(&self) -> String {
        let mut output = String::new();
        let mut begin = 0;
        while begin < self.operations.len() {
            let operation = self.operations[begin];
            let mut end = begin + 1;
            while end < self.operations.len() && self.operations[end] == operation {
                end += 1;
            }
            output.push_str(&format!(
                "{}{}",
                end - begin,
                if operation == b'M' {
                    '='
                } else {
                    operation as char
                }
            ));
            begin = end;
        }
        output
    }
}

#[derive(Debug)]
struct Front {
    origin: i32,
    lo: i32,
    hi: i32,
    values: Vec<i32>,
}

impl Front {
    fn new(lo: i32, hi: i32) -> Result<Self, String> {
        let length = hi
            .checked_sub(lo)
            .and_then(|n| n.checked_add(1))
            .ok_or("WFA diagonal range overflow")?;
        if length < 0 {
            return Err("Invalid WFA diagonal range".into());
        }
        Ok(Self {
            origin: lo,
            lo,
            hi,
            values: super::ksw::filled(length as usize, NULL)?,
        })
    }
    fn null(&self) -> bool {
        self.lo > self.hi
    }
    fn get(&self, k: i32) -> i32 {
        if k < self.lo || k > self.hi {
            return NULL;
        }
        // Recovered Source v2.3 resets initialized bounds to lo/hi after
        // heuristic pruning. compute_init_ends fills any requested padding
        // beyond those active bounds with NULL before recurrence loads;
        // backtrace also explicitly checks lo/hi. Do not resurrect retained
        // pre-pruning offsets as newer canonical WFA2 permits.
        k.checked_sub(self.origin)
            .and_then(|at| usize::try_from(at).ok())
            .and_then(|at| self.values.get(at))
            .copied()
            .unwrap_or(NULL)
    }
    fn set(&mut self, k: i32, offset: i32) {
        self.values[(k - self.origin) as usize] = offset;
    }
    fn equate(&mut self, lo: i32, hi: i32) {
        self.lo = self.lo.max(lo);
        self.hi = self.hi.min(hi);
    }
}

#[derive(Debug)]
struct Level {
    m: Front,
    i: Front,
    d: Front,
}

#[derive(Clone, Copy)]
enum Component {
    M,
    I,
    D,
}

fn front(levels: &[Option<Level>], score: i32, component: Component) -> Option<&Front> {
    let level = usize::try_from(score)
        .ok()
        .and_then(|s| levels.get(s))?
        .as_ref()?;
    let front = match component {
        Component::M => &level.m,
        Component::I => &level.i,
        Component::D => &level.d,
    };
    (!front.null()).then_some(front)
}

fn value(levels: &[Option<Level>], score: i32, component: Component, k: i32) -> i32 {
    front(levels, score, component).map_or(NULL, |wavefront| wavefront.get(k))
}

fn bounds(front: Option<&Front>) -> (i32, i32) {
    // Source's shared null wavefront has exactly these diagonal bounds.
    front.map_or((1, -1), |wavefront| (wavefront.lo, wavefront.hi))
}

fn compute(
    levels: &[Option<Level>],
    score: i32,
    mismatch: i32,
    opening: i32,
    extension: i32,
    pattern: i32,
    text: i32,
) -> Result<Option<Level>, String> {
    let mismatch_score = score - mismatch;
    let open_score = score - opening - extension;
    let extension_score = score - extension;
    let m = front(levels, mismatch_score, Component::M);
    let o = front(levels, open_score, Component::M);
    let i = front(levels, extension_score, Component::I);
    let d = front(levels, extension_score, Component::D);
    if m.is_none() && o.is_none() && i.is_none() && d.is_none() {
        return Ok(None);
    }
    let (ml, mh) = bounds(m);
    let (ol, oh) = bounds(o);
    let (il, ih) = bounds(i);
    let (dl, dh) = bounds(d);
    let lo = ml.min(ol - 1).min(il + 1).min(dl - 1);
    let hi = mh.max(oh + 1).max(ih + 1).max(dh - 1);
    let mut output = Level {
        m: Front::new(lo, hi)?,
        i: Front::new(lo, hi)?,
        d: Front::new(lo, hi)?,
    };
    for k in lo..=hi {
        let inserted = value(levels, open_score, Component::M, k - 1).max(value(
            levels,
            extension_score,
            Component::I,
            k - 1,
        )) + 1;
        let deleted = value(levels, open_score, Component::M, k + 1).max(value(
            levels,
            extension_score,
            Component::D,
            k + 1,
        ));
        let mismatched = value(levels, mismatch_score, Component::M, k) + 1;
        let offset = deleted.max(mismatched.max(inserted));
        output.i.set(k, inserted);
        output.d.set(k, deleted);
        let vertical = offset - k;
        output.m.set(
            k,
            if offset < 0 || offset > text || vertical < 0 || vertical > pattern {
                NULL
            } else {
                offset
            },
        );
    }
    Ok(Some(output))
}

#[derive(Default)]
struct Heuristic {
    wait: i32,
    max_sw: Option<i32>,
}

fn distance(front: &Front, k: i32, pattern: i32, text: i32) -> i32 {
    let offset = front.get(k);
    if offset < 0 {
        -NULL
    } else {
        (pattern - (offset - k)).max(text - offset)
    }
}

impl Heuristic {
    fn cut(&mut self, level: &mut Level, score: i32, pattern: i32, text: i32) {
        let m = &mut level.m;
        if m.null() {
            return;
        }
        self.wait -= 1;
        let old_lo = m.lo;
        let old_hi = m.hi;
        // Source sets the same wait counter in both heuristics: WFadaptive
        // can make XDrop skip this step after it adjusts the wavefront.
        if self.wait <= 0 && m.hi - m.lo + 1 >= 10 {
            let mut nearest = pattern.max(text);
            for k in m.lo..=m.hi {
                nearest = nearest.min(distance(m, k, pattern, text));
            }
            let alignment_k = text - pattern;
            let top = alignment_k.min(m.hi);
            while m.lo < top && distance(m, m.lo, pattern, text) - nearest > 50 {
                m.lo += 1;
            }
            let bottom = alignment_k.max(m.lo);
            while m.hi > bottom && distance(m, m.hi, pattern, text) - nearest > 50 {
                m.hi -= 1;
            }
            self.wait = 1;
        }
        if self.wait <= 0 {
            let sw = |front: &Front, k: i32| {
                let h = front.get(k);
                (-(h - k) - h - score) / 2
            };
            let mut current_max = i32::MIN;
            for k in m.lo..=m.hi {
                if m.get(k) >= 0 {
                    current_max = current_max.max(sw(m, k));
                }
            }
            if let Some(maximum) = self.max_sw {
                while m.lo <= m.hi {
                    if m.get(m.lo) >= 0 && maximum - sw(m, m.lo) < 100 {
                        break;
                    }
                    m.lo += 1;
                }
                while m.hi >= m.lo {
                    if m.get(m.hi) >= 0 && maximum - sw(m, m.hi) < 100 {
                        break;
                    }
                    m.hi -= 1;
                }
                self.max_sw = Some(maximum.max(current_max));
            } else {
                self.max_sw = Some(current_max);
            }
            self.wait = 1;
        }
        if m.lo != old_lo || m.hi != old_hi {
            level.i.equate(m.lo, m.hi);
            level.d.equate(m.lo, m.hi);
        }
    }
}

fn candidate(
    levels: &[Option<Level>],
    score: i32,
    component: Component,
    k: i32,
    increment: i32,
    tag: u8,
) -> Option<(i32, u8)> {
    let offset = value(levels, score, component, k);
    (offset >= 0).then_some((offset + increment, tag))
}

fn append(output: &mut Vec<u8>, operation: u8, count: i32) -> Result<(), String> {
    let count = usize::try_from(count).map_err(|_| "Negative WFA traceback match count")?;
    output
        .try_reserve(count)
        .map_err(|error| error.to_string())?;
    output.resize(output.len() + count, operation);
    Ok(())
}

fn backtrace(
    levels: &[Option<Level>],
    mut score: i32,
    mut k: i32,
    mut offset: i32,
    pattern: i32,
    text: i32,
    mismatch: i32,
    opening: i32,
    extension: i32,
) -> Result<Vec<u8>, String> {
    let mut reverse = Vec::new();
    append(&mut reverse, b'D', pattern - (offset - k))?;
    append(&mut reverse, b'I', text - offset)?;
    let mut component = Component::M;
    while offset - k > 0 && offset > 0 && score > 0 {
        let mismatch_score = score - mismatch;
        let open_score = score - opening - extension;
        let extend_score = score - extension;
        let mut candidates = [None; 5];
        if matches!(component, Component::M) {
            candidates[0] = candidate(levels, mismatch_score, Component::M, k, 1, 9);
        }
        if matches!(component, Component::M | Component::I) {
            candidates[1] = candidate(levels, open_score, Component::M, k - 1, 1, 1);
            candidates[2] = candidate(levels, extend_score, Component::I, k - 1, 1, 2);
        }
        if matches!(component, Component::M | Component::D) {
            candidates[3] = candidate(levels, open_score, Component::M, k + 1, 0, 5);
            candidates[4] = candidate(levels, extend_score, Component::D, k + 1, 0, 6);
        }
        // Source packs the offset above a four-bit type tag: mismatch wins
        // equal offsets, then D-extension/D-open, I-extension/I-open.
        let Some((source, tag)) = candidates.into_iter().flatten().max() else {
            break;
        };
        if matches!(component, Component::M) {
            append(&mut reverse, b'M', offset - source)?;
            offset = source;
            if offset - k <= 0 || offset <= 0 {
                break;
            }
        }
        match tag {
            9 => {
                reverse.push(b'X');
                score = mismatch_score;
                offset -= 1;
                component = Component::M;
            }
            1 | 2 => {
                reverse.push(b'I');
                score = if tag == 1 { open_score } else { extend_score };
                k -= 1;
                offset -= 1;
                component = if tag == 1 { Component::M } else { Component::I };
            }
            5 | 6 => {
                reverse.push(b'D');
                score = if tag == 5 { open_score } else { extend_score };
                k += 1;
                component = if tag == 5 { Component::M } else { Component::D };
            }
            _ => return Err("Invalid WFA traceback tag".into()),
        }
    }
    let mut vertical = offset - k;
    let mut horizontal = offset;
    if matches!(component, Component::M) {
        let matches = vertical.min(horizontal).max(0);
        append(&mut reverse, b'M', matches)?;
        vertical -= matches;
        horizontal -= matches;
        append(&mut reverse, b'D', vertical)?;
        append(&mut reverse, b'I', horizontal)?;
    } else if vertical != 0 || horizontal != 0 || score != 0 {
        return Err("WFA gap-component traceback did not reach the origin".into());
    }
    reverse.reverse();
    Ok(reverse)
}

fn trim(operations: &mut Vec<u8>, mismatch: i32, opening: i32, extension: i32) -> bool {
    let mut maximum = 0;
    let mut maximum_at = 0;
    let mut score = 0;
    let mut last = 0;
    for (at, &operation) in operations.iter().enumerate() {
        score += match operation {
            b'M' => 1,
            b'X' => -mismatch,
            b'I' | b'D' => -extension - if operation == last { 0 } else { opening },
            _ => unreachable!(),
        };
        last = operation;
        if maximum < score {
            maximum = score;
            maximum_at = at;
        }
    }
    let trimmed = maximum_at as i64 != operations.len() as i64 - 1;
    if maximum == 0 {
        operations.clear();
    } else {
        operations.truncate(maximum_at + 1);
    }
    trimmed
}

/// Actual gap-affine score-indexed wavefront recurrence used by Source's
/// WFAlignerGapAffine(0,mismatch,opening,extension,Alignment). Inputs are bytes;
/// Source's wrapper supplies target as pattern and query as text.
pub fn align(
    pattern: &[u8],
    text: &[u8],
    mismatch: i32,
    opening: i32,
    extension: i32,
    extension_mode: bool,
) -> Result<Alignment, String> {
    if mismatch <= 0 || opening < 0 || extension <= 0 {
        return Err(format!(
            "WFA penalties (X={mismatch},O={opening},E={extension}) must be (X>0,O>=0,E>0)"
        ));
    }
    let plen = i32::try_from(pattern.len()).map_err(|_| "WFA pattern exceeds Loc range")?;
    let tlen = i32::try_from(text.len()).map_err(|_| "WFA text exceeds Loc range")?;
    let scope = opening
        .checked_add(extension)
        .ok_or("WFA penalty overflow")?
        .max(mismatch)
        .checked_add(1)
        .ok_or("WFA score scope overflow")?;
    let mut initial = Level {
        m: Front::new(0, 0)?,
        i: Front::new(1, 0)?,
        d: Front::new(1, 0)?,
    };
    initial.m.set(0, 0);
    let mut levels = vec![Some(initial)];
    let mut heuristic = Heuristic {
        wait: 1,
        max_sw: None,
    };
    let mut score = 0i32;
    let mut null_steps = 0;
    loop {
        if let Some(level) = levels[score as usize].as_mut() {
            let wavefront = &mut level.m;
            let mut endpoint = None;
            for k in wavefront.lo..=wavefront.hi {
                let mut offset = wavefront.get(k);
                if offset == NULL {
                    continue;
                }
                let mut vertical = offset - k;
                while offset >= 0
                    && vertical >= 0
                    && offset < tlen
                    && vertical < plen
                    && pattern[vertical as usize] == text[offset as usize]
                {
                    offset += 1;
                    vertical += 1;
                }
                wavefront.set(k, offset);
                let complete = if extension_mode {
                    offset >= tlen || vertical >= plen
                } else {
                    k == tlen - plen && offset >= tlen
                };
                if complete {
                    endpoint = Some((k, offset));
                    break;
                }
            }
            if let Some((k, offset)) = endpoint {
                let mut operations = backtrace(
                    &levels, score, k, offset, plen, tlen, mismatch, opening, extension,
                )?;
                let partial = extension_mode && trim(&mut operations, mismatch, opening, extension);
                return Ok(Alignment {
                    status: i32::from(partial),
                    wavefront_score: score,
                    operations,
                });
            }
            heuristic.cut(level, score, plen, tlen);
        } else if null_steps > scope {
            // With only active WFadaptive/XDrop, Source has no recorded
            // Z-drop endpoint and returns an empty partial CIGAR here.
            return Ok(Alignment {
                status: 1,
                wavefront_score: score,
                operations: Vec::new(),
            });
        }
        score = score.checked_add(1).ok_or("WFA alignment score overflow")?;
        let next = compute(&levels, score, mismatch, opening, extension, plen, tlen)?;
        if next.is_none() {
            null_steps += 1;
        } else {
            null_steps = 0;
        }
        levels.try_reserve(1).map_err(|error| error.to_string())?;
        levels.push(next);
    }
}
