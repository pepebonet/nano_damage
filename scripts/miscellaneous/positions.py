
rounded_outside = [-10, -21, -31, -41, -51, -62,
                   10, 21, 31, 41, 51, 62]

rounded_inside = [-15, -26, -36, -46, -57, -67,
                  15, 26, 36, 46, 57, 67]

ALL_MINORS_IN = []
ALL_MINORS_OUT = []

for r in rounded_inside:
    for i in range(r - 1, r + 2):
        ALL_MINORS_IN.append(i)

for r in rounded_outside:
    for i in range(r - 1, r + 2):
        ALL_MINORS_OUT.append(i)


def interval_minors():
    minors_out = []
    minors_in = []

    for i in range(1, 8):
        current_pos = 10.3 * i
        next_boundary = current_pos + 10.3 / 2
        middle = current_pos + (next_boundary - current_pos) / 2
        minors_out.append([current_pos - (next_boundary - current_pos) / 2, middle])
        minors_in.append([middle, next_boundary + 10.3 / 4])
        minors_out.append([-1 * middle, -1 * (current_pos - \
            (next_boundary - current_pos) / 2), ])
        minors_in.append([-1 * (next_boundary + 10.3 / 4), -1 * middle, ])

    return minors_in, minors_out

MINORS_IN, MINORS_OUT = interval_minors()


DYAD_X_SMALL = [-10.3, -20.6, -30.9, -41.2, -51.5, 10.3, 20.6, 30.9, 41.2, 51.5]
DYAD_X = [-10.3, -20.6, -30.9, -41.2, -51.5, -61.8, -72.1, \
    10.3, 20.6, 30.9, 41.2, 51.5, 61.8, 72.1]
