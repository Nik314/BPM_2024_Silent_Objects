
import pm4py
import copy

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def get_oc_properties(relations):

    convergent_object_types = {a:set() for a in relations["ocel:activity"].unique()}
    divergent_object_types = {a: set() for a in relations["ocel:activity"].unique()}
    deficient_object_types = {a: set() for a in relations["ocel:activity"].unique()}
    related_object_types = {a: set(relations["ocel:type"].unique()) for a in relations["ocel:activity"].unique()}

    look_up_dict_activities = relations.set_index("ocel:eid").to_dict()["ocel:activity"]
    look_up_dict_objects = relations.set_index("ocel:oid").to_dict()["ocel:type"]

    identifiers = relations.groupby("ocel:eid").apply(lambda
        frame:tuple(sorted(set(frame["ocel:oid"].values)))).to_frame(name="all")
    identifiers["activity"] = [look_up_dict_activities[event_id] for event_id in identifiers.index]

    for activity in relations["ocel:activity"].unique():
        sub_relations = relations[relations["ocel:activity"] == activity]
        for object_type in relations["ocel:type"].unique():
            sub_sub_relations = sub_relations[sub_relations["ocel:type"] == object_type]
            if sub_sub_relations["ocel:eid"].nunique() != sub_relations["ocel:eid"].nunique():
                if sub_sub_relations["ocel:eid"].nunique() > 0:
                    deficient_object_types[activity].add(object_type)
                else:
                    related_object_types[activity].remove(object_type)

    for object_type in relations["ocel:type"].unique():
        identifiers[object_type] = identifiers["all"].apply(lambda
            object_set:tuple(sorted(list({object_id for object_id in object_set if look_up_dict_objects[object_id] == object_type}))))

    for object_type in relations["ocel:type"].unique():
        sub_identifiers = identifiers[identifiers[object_type] != set()]
        for activity in relations["ocel:activity"].unique():
            sub_sub_identifiers = sub_identifiers[sub_identifiers["activity"] == activity]

            matches = sub_sub_identifiers.groupby(object_type).apply(lambda
                    frame: frame["all"].nunique())
            matches = matches[[index for index in matches.index if index]]

            if sub_sub_identifiers[object_type].apply(lambda object_set: len(object_set)).max() > 1:
                convergent_object_types[activity].add(object_type)
            if matches.max() > 1:
                divergent_object_types[activity].add(object_type)

    return (convergent_object_types, divergent_object_types,
            deficient_object_types, related_object_types)



def is_combi_identifier(relations, key_set, p, true_identifiers):
    relations = relations[relations["ocel:type"].isin(list(key_set))]
    keys_identifiers = relations.groupby("ocel:eid").apply(lambda frame: tuple(sorted(frame["ocel:oid"]))).nunique()
    return (keys_identifiers / true_identifiers) >= p


def get_minimal_identifiers(relations, p):
    k_list = [k for k in range(1, relations["ocel:type"].nunique() + 1)]
    key_sets = {(t,) for t in relations["ocel:type"].unique()}
    true_identifiers = relations.groupby("ocel:eid").apply(lambda frame: tuple(sorted(frame["ocel:oid"]))).nunique()
    result = set()

    for k in k_list:
        for key_set in copy.deepcopy(key_sets):
            if len(key_set) == k:
                if is_combi_identifier(relations, key_set, p, true_identifiers):
                    result.add(tuple(sorted(key_set)))
        for key_set in copy.deepcopy(key_sets):
            if len(key_set) == k:
                if not key_set in result:
                    for key in [t for t in relations["ocel:type"].unique() if t not in key_set]:
                        new_key_set = set(key_set)
                        new_key_set.add(key)
                        if not any(all(t in new_key_set for t in r) for r in result):
                            key_sets.add(tuple(sorted(new_key_set)))

    return [key for key in result if
            all(key == other_key or len(set(other_key) & set(key)) != len(set(other_key)) for other_key in result)]


def check_shared_identifier(relations, activities, key_set, q):
    relations = relations[relations["ocel:activity"].isin(activities)]
    sub_relations = relations[relations["ocel:type"].isin(key_set)]
    identifiers = [set(sub_relations[sub_relations["ocel:activity"] == a].groupby("ocel:eid").apply(
        lambda frame: tuple(sorted(frame["ocel:oid"]))).values) for a in activities]
    return (len(identifiers[0] & identifiers[1]) / len(identifiers[0] | identifiers[1])) >= q


def get_shared_keys(relations, p, q):
    activity_keys, result = {}, []
    alphabet = list(relations["ocel:activity"].unique())
    for activity in alphabet:
        activity_keys[activity] = get_minimal_identifiers(relations[relations["ocel:activity"] == activity], p)

    for keyset in set(sum(activity_keys.values(), [])):
        relevent_activities = [a for a in alphabet if keyset in activity_keys[a]]
        matrix = {(a, b): 0 for a in relevent_activities for b in relevent_activities}
        done = []
        for a in relevent_activities:
            for b in relevent_activities:
                if (a, b) in done or (b, a) in done:
                    continue
                check = int(check_shared_identifier(relations, [a, b], keyset, q))
                matrix[(a, b)] = check
                matrix[(b, a)] = check
                done.append((a, b))

        graph = csr_matrix([[matrix[(a, b)] for a in relevent_activities] for b in relevent_activities])
        n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
        for i in range(0, n_components):
            result.append({"A": tuple(
                set([relevent_activities[j] for j in range(0, len(relevent_activities)) if labels[j] == i])),
                           "C": keyset})
    return result


def get_invisible_object_types(relations, p, q):

    key_list = get_shared_keys(relations, p, q)
    con,div,dif,rel = get_oc_properties(relations)
    valid = {a:(((rel[a] - div[a]) - con[a]) - dif[a]) for a in rel.keys()}
    result = [key for key in key_list if
              not any([all(ot in valid[a] for a in key["A"]) and
                       set(relations[relations["ocel:type"] == ot]["ocel:activity"].unique()) == set(
                  a for a in key["A"]) for ot in relations["ocel:type"].unique()])]
    return len(result)



ocel = pm4py.read_ocel("your log path.jsonocel")
import time
print(ocel.relations["ocel:oid"].nunique())
print(ocel.relations["ocel:type"].nunique())
print(ocel.relations["ocel:eid"].nunique())
print(ocel.relations["ocel:activity"].nunique())
start = time.time()
result = get_invisible_object_types(ocel.relations,1,1)
print(time.time()-start)
print(result)
