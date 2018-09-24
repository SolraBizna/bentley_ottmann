This is a Rust implementation of the [Bentley­Ottmann algorithm] for finding all intersection points within a set of line segments.

[Bentley­Ottmann algorithm]: https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm

This implementation currently uses 64-bit integers for endpoint coordinates, and ratios of two 64-bit integers for intersection coordinates. These were chosen as adequate for use with FreeType outlines. A future version of this crate will use generics to allow arbitrary suitable types for endpoint coordinates and intersection coordinates.

This implementation is vulnerable to overflow with very large input coordinates. Most likely, even the full range of a 32-bit integer is too much. The original application this was intended for—font rendering—falls well within these limits.

# License

This software is licensed under the zlib license. See [`LICENSE.md`](LICENSE.md) for more information; the short version is that you can do anything you want except claim you wrote this specific software.
