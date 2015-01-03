def motzkinlengthk(length, width, memoize={0: 1}):
    """
    This program computes the number
    of possible Motzkin paths of a specified 
    length and minimum width.

    Motzkin paths can be expressed as a function
    of the Catalan numbers. Catalan numbers are
    represented by sets of parenthesis.
    There are two rules for Catalan numbers:

    1) At any point in the sequence, there
    must be at least as many open parenthesis
    as closed parenthesis

    2) At the end of the sequence, all
    parentheses have been closed (that is,
    there are an equal number of open and
    closed parenthesis).

    This essentialy is the same as the rules for
    using parenthesis in mathematics:

    All parentheses must be paired, and we
    can't have a closed parenthesis that aren't
    preceeded by an open one.

    A Motzkin path is like a Catalan number,
    exept that it contains both parenthesis
    and dots. The parenthesis are all matched
    up, in the same way. We have paths like:
    '(...)(())...'

    The minimum width of a Motzkin path
    is equal to the minimum number of dots
    between any pair of parenthesis.

    This Motzkin path has minimum width 0:
    '(...)(())...'

    This one has minimum width 5:
    '.(.(.....))'

    Motzkin paths can be used represent
    the way in which RNA folds. When RNA folds,
    some of the bases are paired and some
    are not. the paired bases are represented
    by parenthesis. The unpaired bases are
    represented by dots.

    The 'turning/bending radius' of an RNA
    strand is 3 base pairs. That means that
    between any matched base pairs, there
    must be at least three dots. This corresponds
    to a width of 3.

    So, ".(.(...))" works, however,

    ".(.(..))"
    is a violation: The minimum with of this sequence
    is 2, because there are only two dots between the
    innermost parenthesis.

    ".)(.()"
    is a violation, because the first closed parenthesis
    comes before the first open parenthesis


    and ".(.(...)"
    is a violation, as there are more open
    parenthess than closed parenthesis.

    Motzkin paths of minimum length 3 are
    described in the Online Encyclopedia of
    integer sequences here
    https://oeis.org/A023421
    """

    if length not in memoize:
        sequence1shorter = motzkinlengthk(length-1, width)
        sequencecount = 0
        for length_inside_first_set in range(width, (length-2)+1):
            stuff_in_first_set = motzkinlengthk(length_inside_first_set,
                                                width)
            stuff_after_first_set = motzkinlengthk(
                length-2-length_inside_first_set, width)
            sequencecount = sequencecount + \
                stuff_in_first_set * \
                stuff_after_first_set
        memoize[length] = sequence1shorter + sequencecount
    return memoize[length]


# ------------------------------

def motzkin_sequences(length, minwidth):
    """
    This generates all possible Motzkin paths of length
    'length' and minimum width of 'minwidth'

    Call it like this:
    >>> import random
    >>> my_motzkin = motzkin_sequences(8,3)
    >>> all_sequences = [sequence for sequence in my_motzkin]
    >>> len(all_sequences)
    16
    >>> len(set(all_sequences))
    16
    >>> '(.(...))' in all_sequences
    True
    >>> '()......' in all_sequences
    False
    >>> randlen = random.randint(3,15)
    >>> randwidth = random.randint(0, randlen-2)
    >>> my_motzkin = motzkin_sequences(randlen, randwidth)
    >>> genlength = 0
    >>> for i in my_motzkin: genlength += 1
    >>> genlength == motzkinlengthk(randlen, randwidth)
    True
    """
    if length == 0:
        yield ''
    else:
        # Either the first element of this one is a dot, or it isn't
        # If it is a dot, the number of possible sequences is equal
        # to the number of sequences of length 'length-1'
        # (that is, each of these sequences preceeded by a dot)
        dot_before_prev_possibilities = ('.' + seq
                                         for seq
                                         in motzkin_sequences(length-1,
                                                              minwidth))
        for sequence in dot_before_prev_possibilities:
            yield sequence

        # If the first element is not a dot, it is an open parenthesis

        for length_inside_first_set in range(minwidth, length-2+1):
            # If the first element is an open parenthesis, we need to
            # place a matching, closed parenthesis after it.
            # There must be at least minwidth elements between
            # the first parenthesis and its partner, because we
            # need at least minwidth dots between these two parenthesis.
            #
            # There can be as many as length-2 spots between the
            # first and last element, since this would be representative
            # of having the first space be an open parenthesis and the last
            # spot being its partnered close parenthesis.
            #
            # So, let's allocate the parenthesis inside this pair
            length_after_first_set = length-2-length_inside_first_set
            stuff_inside_first_set = [sequence for sequence in
                                      motzkin_sequences
                                      (length_inside_first_set,
                                       minwidth)]
            stuff_after_first_set = [sequence for sequence in
                                     motzkin_sequences(
                                         length_after_first_set,
                                         minwidth)]
            for inside in stuff_inside_first_set:
                for after in stuff_after_first_set:
                    sequence = '({inside}){after}'.format(
                        inside=inside,
                        after=after)
                    yield sequence

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
