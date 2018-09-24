/// A priority queue. We can't use the built-in `BinaryHeap` because we need to
/// be able to remove entries from the queue, and we can't use `BTreeMap` and
/// relatives directly both because they're clumsy for this purpose and because
/// they don't support multiple entries with the same key.

use std::collections::BTreeMap;

pub struct PriorityQueue<K: Copy + Ord, V> {
    inner: BTreeMap<K, Vec<V>>
}

impl<K: Copy + Ord, V> PriorityQueue<K, V> {
    pub fn new() -> PriorityQueue<K, V> {
        PriorityQueue { inner: BTreeMap::new() }
    }
    pub fn insert(&mut self, key: K, value: V) {
        let vec = self.inner.entry(key).or_insert_with(Vec::new);
        vec.push(value);
    }
    pub fn pop(&mut self) -> Option<(K,V)> {
        let key = if let Some(key) = self.inner.keys().next() { *key }
        else { return None };
        let (val, inner_remove) = {
            let vec = self.inner.get_mut(&key).unwrap();
            let val = vec.pop().unwrap();
            (val, vec.is_empty())
        };
        if inner_remove {
            self.inner.remove(&key);
        }
        Some((key, val))
    }
    pub fn remove<H: Fn(&V)->bool>(&mut self, key: &K, checker: H) {
        let inner_remove = if let Some(vec) = self.inner.get_mut(&key) {
            let mut n = 0;
            while n < vec.len() {
                if checker(&vec[n]) {
                    vec.remove(n);
                }
                else { n += 1 }
            }
            vec.is_empty()
        } else { false };
        if inner_remove {
            self.inner.remove(&key);
        }
    }
}
