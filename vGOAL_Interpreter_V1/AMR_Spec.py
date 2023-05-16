import Interpreter as DG
import time

start = time.time()

def main():
    knowledge_base =["forall w. on(w,4) implies available(w)",
                     "forall w. on(w,3) implies available(w)",
                     "forall w. exists p. on(w,p) and equal(p,2) implies delivered(p,w)",
                     "exists p. at(p) and equal(p,6) implies located(charging)",
                     "exists p. at(p) and equal(p,7) implies located(charging)",
                     "exists p. at(p) and equal(p,8) implies located(charging)",
                     "exists l. battery(l) and equal(l,1) implies safe1",
                     "exists l. battery(l) and equal(l,2) implies safe1",
                     "exists p. at(p) and not at(9) implies safe2",
                     #Error implication
                     "E1 implies nonfatal",
                     "E2 implies nonfatal",
                     "E3 implies nonfatal",
                     "E4 implies fatal",
                     "equal(1,1)",
                     "equal(2,2)",
                     "equal(3,3)",
                     "equal(4,4)",
                     "equal(5,5)",
                     "equal(6,6)",
                     "equal(7,7)",
                     "equal(8,8)",
                     "equal(9,9)",
                     "equal(10,10)",
                     "equal(A1,A1)",
                     "equal(A2,A2)",
                     "equal(A3,A3)",
                     "equal(_,_)",
                     "equal(charging,charging)"]
    constraints_of_action_generation = [
        #Ensure the decision-making module will not generate decisions before the revision of goals and beliefs.
        "forall p. at(p) and fatal implies M(p)",
        "nonfatal implies Dummy",
        "forall w,y in D2 . a-goal holding(w) and docked(p) and not holding(y) and docked(4) and available(w) implies A(w)",
        "forall w,y in D2 . a-goal holding(w) and docked(p) and not holding(y) and docked(3) and available(w) implies A(w)",
        "forall p. a-goal docked(p) and at(p) and equal(p,2) implies H(p)",
        "forall p. a-goal docked(p) and at(p) and equal(p,3) implies H(p)",
        "forall p. a-goal docked(p) and at(p) and equal(p,4) implies H(p)",
        "forall p. a-goal docked(p) and at(p) and equal(p,6) implies H(p)",
        "forall p. a-goal docked(p) and at(p) and equal(p,7) implies H(p)",
        "forall p. a-goal docked(p) and at(p) and equal(p,8) implies H(p)",
        #From P3 to P2
        "exists p,x,w. a-goal at(p) and equal(p,2) and docked(x) and equal(x,4) and holding(w) and assigned(p) implies B(x)",
        # From P4 to P2
        "exists p,x,w. a-goal at(p) and equal(p,2) and docked(x) and equal(x,3) and holding(w) assigned(p) implies B(x)",
        "exists p,x,w. a-goal at(p) and equal(p,2) and docked(x) and holding(w) implies S(p)",
        #AMR goes from P1 to P3 or P4, from P6, P7, P8 to P3 or P4, from P2 to P5, from P5 to P6, P7,P8
        "forall p,w in D2 . a-goal at(p) and not holding(w) and not equal(p,2) and assigned(p) implies C(p)",
        "forall p,w in D2 . a-goal at(p) and not holding(w) and not equal(p,2) implies S(p)",
        "forall w. a-goal on(w,2) and at(2) and docked(2) and holding(w) implies D(w,2)",
        "exists l,p. a-goal battery(l) and equal(2,l) and assigned(p) and battery(1) and docked(p) implies E(p)",
        "exists x,y. a-goal located(x) and equal(x,charging) and at(5) implies F",
        "exists x,y. reserved(x,y) implies G"
    ]
    enableness_of_actions = [
        "forall w. A(w) implies pickup(w)",
        "forall p. H(p) implies dock(p)",
        #from P3,P4 to P2
        "exists x. B(x) implies move3(x,2)",
        #from P1 to P2,P3,P4,P6,P7,P8
        "forall p. exists y. C(p) and at(y) and equal(y,1) and not equal(p,5) implies move1(y,p)",
        #from P2,P3,P4,P6,P7,P8 to P5
        "forall p. exists y. C(p) and at(y) and docked(y) and not equal(y,1) and equal(p,5) implies move4(y,p)",
        # from P2,P3,P4,P6,P7,P8 to P5 Docking Error handling
        "forall p. exists y. C(p) and at(y) and not equal(y,1) and equal(p,5) implies move5(y,p)",
        # from P2,P3,P4,P6,P7,P8 to P2,P3,P4,P6,P7,P8
        "forall p. exists y. C(p) and at(y) and not equal(y,1) and not equal(y,5) and not equal(p,1) and not equal(p,5) implies move3(y,p)",
        #from P5 to P6,P7,P8
        "forall p. exists y. C(p) and at(y) and equal(y,5) and implies move2(y,p)",
        "forall w. D(w,2) implies dropoff(w,2)",
        "forall p. E(p) implies charging(p)"
    ]
    sent_message_update = [
        "F implies send!(C) idle(_)",
        "G implies send?(allother) released(_)",
        "forall p. S(p) implies send!(C) idle(p)",
        "forall y. M(y) implies send!(C) at(y)"
                           ]

    event_processing = [
        #Error handling
        "fatal implies drop all",
        "fatal implies delete all",
        "nonfatal and not goal_change implies drop all",
        "nonfatal and not goal_change implies adopt located(charging)",
        "nonfatal and not goal_change implies adopt at(5)",
        "nonfatal and not goal_change implies insert goal_change",
        "nonfatal and E1 implies delete E1",
        "nonfatal and E2 implies delete E2",
        "nonfatal and E3 implies delete E3",
        #Normal event processing
        "forall z. exists x,y. sent!(x) at(y) and reserved(x,z) and not equal(z,y) implies insert idle(z)",
        "forall x,z. exists y. sent!(x) at(y) and reserved(x,z) and not equal(z,y) implies delete reserved(x,z)",
        "forall x,z in D6 . exists y. sent!(x) idle(y) and idle(y) and not reserved(z,y) implies insert reserved(x,y)",
        "forall x. exists y. sent!(x) idle(y) and reserved(x,y) implies send:(x) assigned(y)",
        "exists x,y. sent!(x) idle(y) and reserved(z,y) and equal(x,z) implies delete idle(y)",
        "forall x,y. sent?(x) released(_) and released(y) implies send:(x) idle(y)",
        "forall y. exists x. sent?(x) released(_) and released(y) implies delete released(y)",
        "forall y. exists x. sent:(x) idle(y) implies insert idle(y)",
        "forall y. exists x. sent:(x) idle(y) implies delete reserved(x,y)",
        "forall x,y,z in D6 ,m in D4 . sent!(x) idle(_) and idle(y) and not reserved(z,y) and not reserved(x,m) and equal(y,6) implies insert reserved(x,y)",
        "forall x,y. sent!(x) idle(_) and idle(y) and reserved(x,y) and equal(y,6) implies send:(x) assigned(y)",
        "forall y. exists x. sent!(x) idle(_) and reserved(x,y) and equal(y,6) implies delete idle(y)",
        "forall x,y,z in D6 ,m in D4 . sent!(x) idle(_) and idle(y) and not reserved(z,y) and not reserved(x,m) and equal(y,7) implies insert reserved(x,y)",
        "forall x,y. sent!(x) idle(_) and idle(y) and reserved(x,y) and equal(y,7) implies send:(x) assigned(y)",
        "forall y. exists x. sent!(x) idle(_) and reserved(x,y) and equal(y,7) implies delete idle(y)",
        "forall x,y,z in D6 ,m in D4 . sent!(x) idle(_) and idle(y) and not reserved(z,y) and not reserved(x,m) and equal(y,8) implies insert reserved(x,y)",
        "forall x,y. sent!(x) idle(_) and idle(y) and reserved(x,y) and equal(y,8) implies send:(x) assigned(y)",
        "forall y. exists x. sent!(x) idle(_) and reserved(x,y) and equal(y,8) implies delete idle(y)",
        "forall y. exists x. sent:(x) assigned(y) implies insert assigned(y)",
        "forall p,w. a-goal delivered(p,w) implies adopt on(w,p)",
        "exists p,w. a-goal delivered(p,w) implies adopt located(charging)",
        "exists x,w,p. a-goal on(w,2) and on(w,p) and at(x) implies adopt at(p)",
        "exists x,w,p. a-goal on(w,2) and on(w,p) and at(x) and equal(x,5) implies adopt at(p)",
        "exists x,w,p. a-goal on(w,2) and on(w,p) and docked(x) and equal(x,6) implies adopt at(p)",
        "exists x,w,p. a-goal on(w,2) and on(w,p) and docked(x) and equal(x,7) implies adopt at(p)",
        "exists x,w,p. a-goal on(w,2) and on(w,p) and docked(x) and equal(x,8) implies adopt at(p)",
        "exists w. a-goal on(w,2) and not at(2) implies adopt at(2)",
        "forall y in D2 . exists w. a-goal on(w,2) and not holding(y) implies adopt holding(w)",
        "forall y in D2 . exists p. a-goal located(x) and equal(x,charging) and docked(2) and not holding(y) implies adopt at(5)",
        "forall p. exists x. a-goal located(x) and equal(x,charging) and at(5) and assigned(p) implies adopt at(p)",
        "exists p. a-goal at(p) and equal(p,6) and battery(1) and at(5) implies adopt battery(2)",
        "exists p. a-goal at(p) and equal(p,7) and battery(1) and at(5) implies adopt battery(2)",
        "exists p. a-goal at(p) and equal(p,8) and battery(1) and at(5) implies adopt battery(2)"
    ]

    action_specification = {
        "pickup": "forall w,p,y in D2 . pickup(w) and not holding(y) and on(w,p) implies holding(w) and not on(w,p)",
        #move1: from P1 to P2, P3, P4, P6, P7, P8
        "move1": "forall x,y. move1(x,y) and at(x) implies at(y) and and docked(y) not at(x)",
        #move2: from P5 to P2, P3, P4, P5, P6, P7, P8
        "move2": "forall x,y. move2(x,y) and at(x) implies at(y) and not at(x) and docked(y) and not assigned(x)",
        #move3: from P2, P3, P4, P5, P6, P7, P8 to P2, P3, P4, P5, P6, P7, P8
        "move3":"forall x,y. move3(x,y) and at(x) and docked(x) implies at(y) and not at(x) and docked(y) and not docked(x) and not assigned(x)",
        #move4: from P2, P3, P4, P5, P6, P7, P8 to P5
        "move4": "forall x,y. move4(x,y) and at(x) and docked(x) implies at(y) and not at(x) and not docked(x) and not assigned(x)",
        # move5: from undocked P2, P3, P4, P5, P6, P7, P8  to P5
        "move5": "forall x,y. move5(x,y) and at(x) implies at(y) and not at(x) and not assigned(x)",
        "dropoff": "forall w. dropoff(w) and holding(w) implies on(w,2) and not holding(w)",
        "charging": "exists p. charging(p) and battery(1) implies battery(2) and not battery(1)"
    }
    domain = { "D2": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
              "D4":["6", "7", "8"],"D6": ["A1", "A2","A3"]}
    constants = ["0", "1", "2","3", "4", "5", "6", "7", "8",  "9", "10","charging","allother", "all", '_',"A1","A2","A3","C","D"]
    belief_base4 = ["idle(2)", "idle(3)", "idle(4)","idle(5)", "reserved(A1,6)", "reserved(A2,7)", "reserved(A3,8)"]
    goal_base1 = ['delivered(2,1)']
    goal_base2 = ["delivered(2,2)"]
    goal_base3 = ['delivered(2,3)']
    goal_base4 = ["delivered(2,4)"]
    goals1=[goal_base1]
    goals2 = [goal_base2]
    goals3 = [goal_base3,goal_base4]
    goals4 = []
    dummy_agents=["C"]
    safety = {"A1": ["safe1","safe2"], "A2": ["safe1","safe2"], "A3": ["safe1","safe2"]}
    A1 = DG.Agent("A1", [], goals1)
    A2 = DG.Agent("A2", [], goals2)
    A3 = DG.Agent("A3", [], goals3)
    C = DG.Agent("C", belief_base4, goals4)
    Agents = [A1,A2,A3,C]
    prior_beliefs=["on(1,3)","on(2,4)","on(3,3)","on(4,3)"]
    agent_test = DG.interpreter(Agents, knowledge_base, constraints_of_action_generation,
                                 enableness_of_actions, action_specification, sent_message_update,
                                 event_processing, domain, constants, dummy_agents, safety,prior_beliefs)

    end = time.time()
    f = open("Record.txt", "w+")
    f.write("The duration time is :" + str(end - start))
    f.close()

main()
