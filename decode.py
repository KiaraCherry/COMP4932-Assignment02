import math

QUANT_TABLE = [
    [16, 11, 10, 16, 24, 40, 51, 61],
    [12, 12, 14, 19, 26, 58, 60, 55],
    [14, 13, 16, 24, 40, 57, 69, 56],
    [14, 17, 22, 29, 51, 87, 80, 62],
    [18, 22, 37, 56, 68, 109, 103, 77],
    [24, 35, 55, 64, 81, 104, 113, 92],
    [49, 64, 78, 87, 103, 121, 120, 101],
    [72, 92, 95, 98, 112, 100, 103, 99]
]

ZIGZAG_COORDS = [
    (0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2), (0, 3), (1, 2),
    (2, 1), (3, 0), (4, 0), (3, 1), (2, 2), (1, 3), (0, 4), (0, 5),
    (1, 4), (2, 3), (3, 2), (4, 1), (5, 0), (6, 0), (5, 1), (4, 2),
    (3, 3), (2, 4), (1, 5), (0, 6), (0, 7), (1, 6), (2, 5), (3, 4),
    (4, 3), (5, 2), (6, 1), (7, 0), (7, 1), (6, 2), (5, 3), (4, 4),
    (3, 5), (2, 6), (1, 7), (2, 7), (3, 6), (4, 5), (5, 4), (6, 3),
    (7, 2), (7, 3), (6, 4), (5, 5), (4, 6), (3, 7), (4, 7), (5, 6),
    (6, 5), (7, 4), (7, 5), (6, 6), (5, 7), (6, 7), (7, 6), (7, 7)
]

RLC_INPUT = [
    -26,
    (0, -3), (1, -3), (0, -2), (0, -6), (0,  2), (0, -4),
    (0,  1), (0, -3), (0,  1), (0,  1), (0,  5), (0,  1),
    (0,  2), (0, -1), (0,  1), (0, -1), (0,  2), (5, -1), (0, -1), (0,  0)
]


def decode_rlc(rlc_input):
    """
    Decode a run-length code sequence into a 64-element 1D vector.

    :param rlc_input: A list of run-length code sequences.
    :return: A 64-element 1D vector.
    """
    vector = [0] * 64
    vector[0] = rlc_input[0]
    pos = 1
    for (run, val) in rlc_input[1:]:
        if run == 0 and val == 0:
            break
        pos += run
        vector[pos] = val
        pos += 1
    return vector


def inverse_zigzag(vector):
    """
    Reconstruct 8x8 block from 64-element vector using standard JPEG zigzag scan order.

    :param vector: A 64-element 1D vector.
    :return: A 8x8 nested list of integers.
    """
    blocks = [[0] * 8 for _ in range(8)]
    for idx, (row, col) in enumerate(ZIGZAG_COORDS):
        blocks[row][col] = vector[idx]
    return blocks


def dequantize(blocks):
    """
    Dequantize a 8x8 block by multiplying with standard JPEG luminance quantization.

    :param blocks: 8x8 nested list.
    :return: 8x8 nested list of dequantized DCT coefficients.
    """
    result = [[0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            result[i][j] = blocks[i][j] * QUANT_TABLE[i][j]
    return result


def c(xi):
    """
    Compute scaling factor.

    :param xi: A non-negative integer.
    :return: The scaling factor as a float.
    """
    return math.sqrt(2) / 2 if xi == 0 else 1.0


def idct_2d(F):
    """
    Perform 2D Inverse Discrete Cosine Transform (IDCT) on an 8x8 block.

    :param F: a 8x8 nested list of dequantized coefficients.
    :return: A 8x8 nested list of integers.
    """
    f = [[0.0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            total = 0.0
            for u in range(8):
                for v in range(8):
                    total += (
                        c(u) * c(v) / 4.0
                        * math.cos((2 * i + 1) * u * math.pi / 16)
                        * math.cos((2 * j + 1) * v * math.pi / 16)
                        * F[u][v]
                    )
            f[i][j] = round(total)
    return f


def level_shift(blocks):
    """
    Apply a level shift by adding 128 to each element in the block.

    :param blocks: a 8x8 nested list of integers.
    :return: A 8x8 nested list of integers.
    """
    return [[blocks[i][j] + 128 for j in range(8)] for i in range(8)]


def print_matrix(matrix, title):
    """
    Print matrix.

    :param matrix: The matrix to print.
    :param title: The title of the matrix.
    """
    print(f"\n{title}")
    print("-" * 50)
    for row in matrix:
        print("  " + "  ".join(f"{val:4d}" for val in row))


def main():
    """
    Drive the program.
    """
    vector = decode_rlc(RLC_INPUT)
    print("Step 1 — Decoded 64-vector:")
    print(" ", vector)

    # Step 2: Inverse Zigzag
    blocks = inverse_zigzag(vector)
    print_matrix(blocks, "Step 2 — 8x8 block after inverse zigzag:")

    # Step 3: Dequantization
    dequant = dequantize(blocks)

    print_matrix(dequant, "Step 3 — 8x8 block after dequantization:")

    # Step 4: 2D IDCT
    idct_result = idct_2d(dequant)
    print_matrix(idct_result, "Step 4 — 8x8 block after IDCT:")

    # Step 5: Level shift (+128)
    final = level_shift(idct_result)
    print_matrix(final, "Step 5 — Final reconstructed pixel values (+128):")

    # Original 8x8 image block
    original = [
        [52, 55, 61, 66, 70, 61, 64, 73],
        [63, 59, 55, 90, 109, 85, 69, 72],
        [62, 59, 68, 113, 144, 104, 66, 73],
        [63, 58, 71, 122, 154, 106, 70, 69],
        [67, 61, 68, 104, 126, 88, 68, 70],
        [79, 65, 60, 70, 77, 68, 58, 75],
        [85, 71, 64, 59, 55, 61, 65, 83],
        [87, 79, 69, 68, 65, 76, 78, 94],
    ]
    print_matrix(original, "Bonus — Original image block (for comparison):")

    # Display images visually using matplotlib
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        fig = plt.figure(figsize=(10, 5))
        fig.suptitle("JPEG Decoder: Reconstructed vs Original", fontsize=14)
        gs = gridspec.GridSpec(1, 2, wspace=0.4)

        ax1 = fig.add_subplot(gs[0])
        ax1.imshow(final, cmap='gray', vmin=0, vmax=255, interpolation='nearest')
        ax1.set_title("Reconstructed")
        ax1.axis('off')

        ax2 = fig.add_subplot(gs[1])
        ax2.imshow(original, cmap='gray', vmin=0, vmax=255, interpolation='nearest')
        ax2.set_title("Original")
        ax2.axis('off')

        plt.savefig("jpegComparison.png", dpi=150, bbox_inches='tight')
        print("\nImage saved as jpegComparison.png")
        plt.show()

    except ImportError:
        print("\nInstall matplotlib to display the visual comparison.")


if __name__ == '__main__':
    main()
